"""
This implements the Cado-NFS api server. Some useful
documentation sources that I perused in order to code this.

https://blog.miguelgrinberg.com/post/running-your-flask-application-over-https
https://stackoverflow.com/questions/39853643/is-it-possible-to-mark-a-method-with-a-route-in-flask

The cado-nfs api server is built upon:

 - flask (python3-flask) in order to actually build an api
 - optionally, gunicorn (python3-gunicorn) in order to make a
   multithreaded server from it. (still WIP)

It serves two rather different audiences.

 - The endpoints that cado-nfs-client.py needs (/workunit, /upload,
   /file, /files, /WUstatus) are unauthenticated, and gated by
   server.whitelist alone. That is the pre-existing design: a client is
   given a url and a certificate fingerprint, never a secret.

 - Everything under /api/v1, plus the web ui at /ui/ that consumes it,
   is meant for whoever runs the computation. Those endpoints require a
   bearer token which we write to <workdir>/<name>.api-token with mode
   0600, and are gated by server.ui_whitelist (localhost by default, so
   that an ssh tunnel works without widening the client whitelist).

The OpenAPI 3.1 description of the whole thing is assembled in
cadofactor/api/spec.py with no third-party tooling, and served at
/api/v1/openapi.json.
"""

import flask
import hashlib
import json
import logging
import multiprocessing
import os.path
import re
import socket
import sys
import threading
import time

from cadofactor.cadofactor_tools import UploadDirProvider
from cadofactor import wudb
from cadofactor.api import auth
from cadofactor.api import spec
from cadofactor.api.spec import api_route
from cadofactor.api.admin import AdminEndpoints, spec_schemas
from cadofactor.api.views import ServingState, DbViews
from cadofactor.api.pool import DbSession, DbSessionPool
from werkzeug import serving

from ipaddress import ip_address, ip_network

import werkzeug
from werkzeug.utils import secure_filename
from werkzeug.exceptions import HTTPException

# See #30142
werkzeug.utils._filename_ascii_strip_re = re.compile(r"[^A-Za-z0-9,_.-]")

try:
    import gunicorn.app.base
except ImportError:
    pass

try:
    # interestingly, it seems that I don't even need the ssl module at all.
    import ssl  # noqa: F401
    from cadofactor.cadofactor_tools.certificate import \
        get_server_alternate_names, \
        create_certificate, \
        get_certificate_hash
except ModuleNotFoundError:
    pass

HAVE_SSL = 'ssl' in sys.modules

# Where the (dependency-free, build-step-free) web ui lives.
UI_DIRECTORY = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'api', 'ui')

API_PREFIX = "/api/v1"

SPEC_INFO = dict(
    title="cado-nfs api",
    version="1.0.0",
    description="Monitoring and administration of a running cado-nfs"
                " computation, plus the workunit endpoints used by"
                " cado-nfs-client.py.",
    license_name="GNU LGPL 2.1",
    license_url="https://www.gnu.org/licenses/old-licenses/"
                "lgpl-2.1.en.html",
    contact={"name": "The Cado-NFS development team",
             "url": "https://cado-nfs.inria.fr"},
)


class ApiServer(flask.Flask):
    """
    This is the cado-nfs api server. Its job is only to handle
    incoming requests, and communicate with the database. In turn, the
    rest of cado-nfs communicates with the database. (the wait loop is in
    cadotask.ClientServerTask.wait)
    """

    def __init__(self,
                 serveraddress,
                 serverport,
                 dbdata,
                 threaded=None,
                 debug=False,
                 uploaddir=None,
                 nrsubdir=None,
                 # scriptdir=None,
                 only_registered=True,
                 cafile=None,
                 whitelist=None,
                 ui_whitelist=None,
                 timeout_hint=None,
                 workdir=None,
                 name=None,
                 read_only=False,
                 # linger_before_quit=False
                 ):
        # some parameters are currently not targeted by the
        # implementation as it is now.

        """
        This starts the server at serveraddress:serverport, using the dbdata
        connection to the database.

        dbdata is a wudb.DBFactory object, which is basically a wrapper
        around an uri. Workers can use to make a fresh
        connection to the database by calling its connect() method. Note
        that dbdata itself is not a connection object!

        The threaded argument can be used to run a production WSGI server
        with gunicorn *if* this software is available. If not, we fall
        back to using a single-thread flask server, which should be
        sufficient for most purposes. Note that threaded=None is used as
        a default behaviour, where a threaded server is used
        opportunistically unless debug mode is on.

        The debug argument spawns a more noisy server.

        There are a few situations where it's desirable for the api server to
        process requests directly without even talking to the db. Namely, we
        want it to be able to return select files from the workdir (factor
        bases and such). Those are pulled from the database.
        The only_registered flag is a priori always equal to True. If False,
        a blanket access to all files is offered, which is probably unwise.

        TLS is used if cafile is given, in which case an ad hoc certificate
        is written to this file name.

        whitelist gates the endpoints that cado-nfs-client.py uses.
        ui_whitelist gates the monitoring and administration endpoints
        under /api/v1 and the web ui at /ui/; it defaults to the loopback
        addresses, so that "ssh -L 8001:localhost:8001" gives access to
        the ui without opening the client endpoints any wider.

        workdir and name locate the files that belong to the
        computation: the api token is written to
        <workdir>/<name>.api-token, and <workdir>/<name>.log is what the
        /api/v1/log endpoint tails. Without them the api server still
        serves clients, but the authenticated part of the api is
        disabled, since there is nowhere to put the token.

        read_only is for looking at a working directory that no
        computation is currently driving (cado-nfs.py --ui-only). The
        monitoring endpoints work as usual, but no workunit is handed
        out and no action is accepted: with no task running, there is
        nobody to pick up a resubmission, and handing leftover
        workunits to a client that happened to connect would be
        actively wrong.
        """

        self.name = "API server"
        super().__init__(self.name)
        self.logger = logging.getLogger(self.name)

        if True:
            # check for #30119
            import importlib
            import importlib.metadata

            def Version(c):
                return tuple([int(x) for x in c.split('.')])

            vw = Version(importlib.metadata.version('werkzeug'))
            v313 = Version("3.1.3")
            if vw <= v313:
                message = r"""
The flask package that is used here is affected by a bug
in the werkzeug package: See:
 - https://github.com/pallets/werkzeug/issues/3065
 - https://gitlab.inria.fr/cado-nfs/cado-nfs/-/issues/30119
It is expected that a future version of werkzeug (more
recent than 3.1.3) fixes this.
                """
                for m in message.strip().split("\n"):
                    self.logger.warning(m)

        self.address = serveraddress if serveraddress else "0.0.0.0"
        self.port = serverport

        self.cafile = cafile

        if not HAVE_SSL:
            if self.cafile:
                self.logger.warning("ssl not available,"
                                    f" cafile={cafile} ignored")
            self.cafile = None
            self._ssl_context = None

        if self.cafile is not None:
            if 'ssl' not in sys.modules:
                raise RuntimeError("python ssl module is missing")
            SAN = get_server_alternate_names(serveraddress)
            self._ssl_context = create_certificate(self.cafile,
                                                   self.address,
                                                   SAN)
            if self._ssl_context is None:
                self.logger.warning("ssl not available,"
                                    f" cafile={cafile} ignored")
                self.cafile = None
        else:
            self._ssl_context = None

        # https://stackoverflow.com/questions/70396641/how-to-run-gunicorn-inside-python-not-as-a-command-line
        if threaded is None:
            threaded = False if debug else True

        # How many requests we intend to have in flight at once. Note
        # that werkzeug's threaded server does not limit its threads to
        # this -- it starts one per request -- so this is a target for
        # sizing the database pool, and the pool is then what actually
        # bounds concurrent database work.
        if not threaded:
            self.nthreads = 1
        elif type(threaded) is int:
            self.nthreads = threaded
        else:
            self.nthreads = min(multiprocessing.cpu_count() * 2, 4) + 1

        # inject fields to the api app object.

        self.database_uri = dbdata
        # One pooled session per server thread, plus a little slack for
        # the odd concurrent request. Bounding it also bounds how many
        # writers queue against sqlite, which is no bad thing.
        self._pool = DbSessionPool(self._new_session,
                                   maxsize=self.nthreads + 2)
        self._offline_session = None
        self._offline_lock = threading.Lock()
        self.wuaccess = None
        self.only_registered = only_registered
        self.upload_dir_provider = UploadDirProvider(uploaddir, nrsubdir)
        self.timeout_hint = timeout_hint
        self.whitelist = self._resolve_whitelist(whitelist)
        self.logger.info(f"server whitelist is {self.whitelist}")
        self.ui_whitelist = self._resolve_whitelist(
            ui_whitelist if ui_whitelist is not None
            else ["127.0.0.1/32", "::1/128"])
        self.logger.info(f"ui whitelist is {self.ui_whitelist}")
        self.wdir = workdir
        self.computation_name = name

        # Whether we are still handing workunits out lives in the
        # database rather than in this object. It is written from
        # cado-nfs.py's main thread and read by request handlers, and
        # those two need not be the same process -- they already are not
        # under gunicorn, whose arbiter forks its workers. Note that
        # get_db_connection is passed as a callable, not called: sqlite3
        # connections are bound to the thread that created them, which
        # is why _get_db_things keeps a per-thread pool.
        self.read_only = read_only
        self.serving = ServingState(self.session)
        # Starting up means we are serving again -- unless we are only
        # here to look, in which case we must not resurrect a flag that
        # an earlier run left cleared.
        self.serving.set(not read_only)

        # Read-only queries behind the monitoring endpoints.
        self.views = DbViews(self.session)

        # The authenticated part of the api needs somewhere private to
        # put its token. Without a workdir we simply do not offer it.
        self.api_token = None
        self.token_file = None
        if workdir is not None and name is not None:
            self.token_file = os.path.join(workdir, name + ".api-token")
            try:
                self.api_token = auth.load_or_create_token(self.token_file)
            except (auth.TokenError, OSError) as e:
                self.logger.error("Could not set up the api token: %s", e)
                self.logger.error("The monitoring api and the web ui will"
                                  " not be available")
                self.token_file = None
        else:
            self.logger.info("No workdir given to the api server;"
                             " the monitoring api and the web ui are"
                             " disabled")

        self.admin = AdminEndpoints(self)
        self._route_endpoints()

        if False and threaded and 'gunicorn' in sys.modules:
            self._prepare_run_gunicorn(serveraddress, serverport, threaded)
        else:
            # self._prepare_run_flask_debug_server(
            #   serveraddress, serverport, threaded, debug=debug)
            self._prepare_run_werkzeug(serveraddress, serverport, threaded)
            log = logging.getLogger('werkzeug')
            log.setLevel(logging.WARNING)

    def _prepare_run_gunicorn(self, serveraddress, serverport, threaded):
        class StandaloneApplication(gunicorn.app.base.BaseApplication):
            def __init__(self, app, options={}):
                self.options = options
                self.application = app
                super().__init__()
                # I'd love to be able to run with gunicorn, but
                # unfortunately I can't find a way to find the port it
                # binds to! The following is from the werkzeug server.
                # self.port = self.server.socket.getsockname()[1]

            def load_config(self):
                config = {key: value
                          for key, value in self.options.items()
                          if key in self.cfg.settings and value is not None
                          }
                for key, value in config.items():
                    self.cfg.set(key.lower(), value)

            def load(self):
                return self.application

        self.logger.info("Running from Gunicorn")

        if type(threaded) is int:
            nthreads = threaded
        else:
            nthreads = min(multiprocessing.cpu_count() * 2, 4) + 1

        def pre_fork(server, worker):
            # for debugging
            print(f"pre-fork server {server} worker {worker}",
                  file=sys.stderr)

        options = {
                'bind': f"{self.address}:{serverport}",
                'workers': nthreads,
                # 'threads': number_of_workers(),
                'timeout': 120,
                'pre_fork': pre_fork,
        }

        if self._ssl_context is not None:
            options['certfile'] = self._ssl_context[0]
            options['keyfile'] = self._ssl_context[1]

        self._run_object = StandaloneApplication(self, options)
        self._run_args = []
        self._run_kwargs = {}

    def _prepare_run_flask_debug_server(self,
                                        serveraddress, serverport,
                                        threaded,
                                        debug=False):
        if threaded:
            self.logger.warning("Note: cannot run threaded server since"
                                " we do not have gunicorn.")
        self.logger.info("Running from flask")
        self._run_object = self
        self._run_args = []
        self._run_kwargs = dict(host=serveraddress,
                                port=serverport,
                                debug=debug,
                                use_reloader=False)
        if self._ssl_context is not None:
            self._run_kwargs['ssl_context'] = self

    def _prepare_run_werkzeug(self,
                              serveraddress, serverport,
                              threaded):
        class Foo(object):
            def __init__(self, address, port, app, **kwargs):
                self.options = options
                self.application = app
                self.server = serving.make_server(host=address,
                                                  port=port,
                                                  app=self.application,
                                                  **kwargs)
                self.port = self.server.socket.getsockname()[1]

            def load_config(self):
                pass

            def load(self):
                return self.application

            def run(self):
                from threading import Thread
                thread = Thread(target=self.server.serve_forever)
                thread.daemon = True
                thread.start()

        options = {}

        nthreads = self.nthreads

        if threaded:
            options['threaded'] = nthreads

        if self._ssl_context is not None:
            options['ssl_context'] = self._ssl_context

        self._run_object = Foo(self.address, serverport, self, **options)
        self._run_args = []
        self._run_kwargs = {}
        self.logger.info("Running from werkzeug (%d thread(s))", nthreads)

        scheme = 'https' if self._ssl_context else 'http'
        bound_address = self._run_object.server.socket.getsockname()[0]
        bound_port = self._run_object.server.socket.getsockname()[1]
        self.url = "%s://%s:%d" % (scheme, bound_address, bound_port)

        self.port = bound_port
        self.url = self.url.replace('0.0.0.0', 'localhost')
        self.server = self._run_object.server

    def serve(self):
        self.logger.info(f"Running on {self.url} (Press CTRL+C to quit))")
        self._run_object.run(*self._run_args, **self._run_kwargs)
        # document the connection procedure.
        connection_parameters = "--server=%s" % self.url
        if self._ssl_context:
            h = get_certificate_hash(self._ssl_context[0])
            connection_parameters += " --certsha1=%s" % h

        self.logger.info("You can start additional cado-nfs-client.py scripts"
                         " with parameters: " + connection_parameters)
        self.logger.info("If you want to start additional clients, remember "
                         "to add their hosts to server.whitelist")

        if self.api_token is not None:
            # The token travels in the fragment, which browsers keep to
            # themselves: it appears in no request line, no access log
            # and no Referer header.
            self.logger.info("Web UI: %s/ui/#token=%s",
                             self.url, self.api_token)
            self.logger.info("  (or open %s/ui/ and paste the token from"
                             " %s)", self.url, self.token_file)
            self.logger.info("  Terminal UI: cado-nfs-monitor.py"
                             " %s --workdir=%s watch",
                             connection_parameters, self.wdir)
            self.logger.info("  The UI listens for %s only; use"
                             " 'ssh -L %d:localhost:%d <host>' from"
                             " elsewhere, or set server.ui_whitelist",
                             ", ".join(self.ui_whitelist) or "nobody",
                             self.port, self.port)

    def get_cert_sha1(self):
        if self._ssl_context:
            return get_certificate_hash(self._ssl_context[0])
        else:
            return None

    def get_url(self, origin=None):
        # if origin=localhost, then we want to insist on using localhost
        # as the peer name. In fairness, I'm not sure it's really
        # important. It's only used in the all-clients-on-localhost case,
        # and then the default listen address of 0.0.0.0 gets rewritten
        # to localhost anyway.
        return self.url

    def get_port(self):
        return self.port

    def stop_serving_wus(self):
        self.logger.info("Got notification to stop serving Workunits")
        self.serving.set(False)

    def resume_serving_wus(self):
        self.logger.info("Resuming serving Workunits")
        self.serving.set(True)

    def shutdown(self, exc=None):
        if exc is not None:
            self.logger.info("Shutting down server on exception %s", exc)
        else:
            self.logger.info("Shutting down server")
        # app.shutdown()
        self.server.shutdown()

    def _new_session(self):
        return DbSession(
            lambda: self.database_uri.connect(shared_across_threads=True))

    def session(self):
        """
        The database session belonging to the work in hand.

        Inside a request it is borrowed from the pool and given back
        when the request ends, so the cost of opening a connection and
        of creating the dictionary tables is paid once per pooled
        session rather than once per request. That matters because
        werkzeug's threaded server runs a *new thread for every
        request*: keyed by thread, as this used to be, every request
        opened a connection and took an EXCLUSIVE lock to create tables
        that already existed, on the very database the computation is
        trying to use.

        Outside a request -- cado-nfs.py's own thread calling
        stop_serving_wus(), say -- a single long-lived session is used
        instead, under a lock.
        """
        if flask.has_request_context():
            session = getattr(flask.g, "cado_session", None)
            if session is None:
                session = self._pool.borrow()
                flask.g.cado_session = session
            return session
        with self._offline_lock:
            if self._offline_session is None:
                self._offline_session = self._new_session()
            return self._offline_session

    def _release_session(self, exc=None):
        session = getattr(flask.g, "cado_session", None)
        if session is not None:
            flask.g.cado_session = None
            # A request that blew up may have left a transaction open;
            # do not hand that to the next borrower.
            self._pool.release(session, discard=exc is not None)

    def get_db_connection(self):
        return self.session().connection

    def get_wuaccess(self):
        return self.session().wuaccess

    def get_registered_filenames(self):
        return self.session().registered_filenames

    def get_upload_folder(self, key):
        """
        Pick the upload subdirectory that an uploaded file goes to.

        The point of the subdirectories is only to keep any single
        directory from growing to a size that the filesystem handles
        badly, so any reasonably uniform spread will do. We derive it
        from the file's own name rather than from a counter, because a
        counter living in this object would spread unevenly as soon as
        the server runs as more than one process.
        """
        h = int.from_bytes(hashlib.sha1(key.encode('utf-8',
                                                   'replace')).digest()[:4],
                           'big')
        return self.upload_dir_provider(h)

    def get_workdir(self):
        """
        Gets the main work directory as it is stored in the database
        """
        d = wudb.DictDbAccess(self.get_db_connection(), 'tasks')
        return d['workdir']

    @staticmethod
    def _resolve_whitelist(entries):
        """
        Turn host names into /32 networks, leave addresses alone.

        Only names are looked up. Handing an address or a CIDR block to
        gethostbyname() asks the resolver a question with no answer,
        and on a machine whose resolver is slow to say so -- a mac with
        search domains and no reachable server, say -- each of those
        non-answers can take seconds. The server used to do three of
        them before it could serve anything at all, which was enough to
        make the test suite time out on such a machine.

        A bare address is left as it is rather than given a /32: it is
        already a one-address network, and api_limit_remote_addr parses
        it with the same ip_network() either way.

        >>> ApiServer._resolve_whitelist(['127.0.0.1/32', '::1/128'])
        ['127.0.0.1/32', '::1/128']
        >>> ApiServer._resolve_whitelist(['10.0.0.1'])
        ['10.0.0.1']
        >>> ApiServer._resolve_whitelist(['0.0.0.0/0'])
        ['0.0.0.0/0']
        >>> ApiServer._resolve_whitelist([])
        []
        """
        resolved = []
        for w in entries or []:
            try:
                # Already an address or a block: nothing to ask anyone.
                ip_network(w, strict=False)
                resolved.append(w)
                continue
            except ValueError:
                pass
            try:
                resolved.append(f"{socket.gethostbyname(w)}/32")
            except OSError:
                resolved.append(w)
        return resolved

    def _route_endpoints(self):
        self.errorhandler(HTTPException)(self.api_errorhandler)
        self.before_request(self.api_limit_remote_addr)
        self.teardown_request(self._release_session)

        # Routes and their documentation come from the same @api_route
        # declaration, so an endpoint cannot be reachable and yet absent
        # from the OpenAPI document. See cadofactor/api/spec.py.
        self.api_routes = []
        for method, meta in (spec.collect_api_routes(self)
                             + spec.collect_api_routes(self.admin)):
            if meta["auth"] and self.api_token is None:
                # No token could be set up, so refuse to expose the
                # endpoints that the token is meant to protect.
                continue
            self.api_routes.append(meta)
            self.route(meta["rule"], methods=list(meta["methods"]))(method)

    def _peer_allowed(self, peer, whitelist):
        try:
            address = ip_address(peer)
        except ValueError:
            return False
        for net in whitelist or []:
            try:
                if address in ip_network(net):
                    return True
            except ValueError:
                self.logger.warning("ignoring malformed whitelist"
                                    " entry %s", net)
        return False

    def is_ui_path(self, path):
        """
        Whether a path belongs to the part of the server that is meant
        for whoever runs the computation, rather than for clients.

        >>> ApiServer.is_ui_path(None, '/api/v1/clients')
        True
        >>> ApiServer.is_ui_path(None, '/ui/app.js')
        True
        >>> ApiServer.is_ui_path(None, '/api/docs')
        True
        >>> ApiServer.is_ui_path(None, '/workunit')
        False
        >>> ApiServer.is_ui_path(None, '/upload')
        False
        """
        return (path.startswith(API_PREFIX)
                or path.startswith("/api/docs")
                or path == "/ui"
                or path.startswith("/ui/"))

    def api_limit_remote_addr(self):
        """
        Implements ip filtering.

        Client endpoints are gated by server.whitelist, as they always
        were. The ui and the monitoring api are gated by
        server.ui_whitelist as well, which defaults to loopback only --
        the point being that "ssh -L" reaches the ui without the client
        whitelist having to be widened to the operator's workstation.
        """
        peer = flask.request.remote_addr
        path = flask.request.path

        if self._peer_allowed(peer, self.whitelist):
            return
        if self.is_ui_path(path) \
                and self._peer_allowed(peer, self.ui_whitelist):
            return

        self.logger.error(f'blocked incoming request from {peer}')
        if not self.whitelist:
            self.logger.error(' NOTE: no whitelist is configured,'
                              ' all ip addresses are blocked anyway.'
                              ' You probably want to add'
                              ' a server.whitelist= argument')
        flask.abort(403)  # Forbidden

    def api_errorhandler(self, e):
        """
        Return JSON instead of HTML for HTTP errors.
        """
        # start with the correct headers and status code from the error
        response = e.get_response()
        # replace the body with JSON
        response.data = json.dumps({
            "code": e.code,
            "name": e.name,
            "description": e.description,
        })
        response.content_type = "application/json"
        # Return the response as it is, keeping the status code that the
        # exception carries. Appending a code here -- as this used to do
        # with a hardcoded 404 -- overrides it, and in particular turned
        # the 410 that /workunit answers at the end of the distributed
        # phase into a 404. cado-nfs-client.py keys its clean shutdown
        # off that 410 (see WorkunitClientToFinish), so it never got the
        # message and kept retrying instead.
        return response

    @api_route("/", tags=["misc"],
               summary="Liveness probe",
               description="Answers as soon as the server is up. Used by"
                           " clients and by monitoring tools to tell a"
                           " running server from a closed port.",
               responses={200: ("The server is alive",
                                {"type": "object",
                                 "properties": {
                                     "message": {"type": "string"}}})})
    def api_hello_world(self):
        resp = {'message': "Hello, World!"}

        return flask.json.jsonify(resp), 200

    @api_route("/workunit", tags=["client"],
               summary="Hand a fresh workunit to a client",
               description="The client must identify itself with the"
                           " clientid form field. A 404 is a priori"
                           " temporary: it only means that the server"
                           " has not provisioned fresh workunits yet.",
               parameters=[
                   spec.query_parameter(
                       "clientid", {"type": "string"},
                       "client-defined identifier", required=True)],
               responses={
                   200: ("A fresh workunit the client can work on",
                         {"$ref": "#/components/schemas/Workunit"}),
                   403: "No clientid was provided",
                   404: "No work available for the time being",
                   410: "The distributed computation phase is over."
                        " Clients should terminate now.",
               })
    def api_get_workunit(self):
        clientid = flask.request.form.get('clientid')
        if clientid is None:
            # self.logger.debug(f"got {flask.request.path}"
            #                   f" from {flask.request.remote_addr}"
            #                   f" (no client identification provided)")
            flask.abort(403, 'clientid must be provided')

        # self.logger.debug(f"got {flask.request.path}"
        #                   f" from {flask.request.remote_addr}"
        #                   f" with client identification {clientid}")

        if self.read_only:
            flask.abort(410, "This server is only serving the monitoring"
                             " interface; no computation is running")
        if not self.serving.get():
            flask.abort(410, "Distributed computation finished")

        # we might want to make the timeout dependent on the clientid.
        wu = self.get_wuaccess().assign(clientid,
                                        timeout_hint=self.timeout_hint)

        if not wu:
            flask.abort(404, "No work available")

        self.logger.info(f"Sending workunit {wu.get_id()}"
                         f" to client {clientid}")

        # change the timeout to a deadline for the json that we return to the
        # client. Inside the database, we really only maintain a timeout
        # value.
        if 'timeout' in wu:
            wu['deadline'] = time.time() + float(wu['timeout'])
            del wu['timeout']

        return flask.json.jsonify(wu), 200

    @api_route("/file/<path:path>", tags=["client"],
               summary="Download a file registered for download",
               description="Clients use this to fetch the input files"
                           " and executables that their workunits refer"
                           " to. Only files that a task has registered"
                           " are served, unless server.only_registered"
                           " is false.",
               responses={200: "The file contents",
                          404: "No such registered file"})
    def api_download_file(self, path):
        d = self.get_registered_filenames()
        if d is None:
            if re.match('^/', path):
                full_path = path
            else:
                full_path = os.path.join(self.get_workdir(), path)
            self.logger.info(f"got request for arbitrary file {path},"
                             f" which resolves to {full_path}")
            if os.path.isfile(full_path):
                return flask.send_file(full_path)
        else:
            file = d.get(path)
            if file is not None:
                self.logger.info(f"got request for file {path},"
                                 f" which resolves to {file}")
                if os.path.isfile(file):
                    dirname, basename = os.path.split(file)
                    return flask.send_from_directory(dirname, basename)
        flask.abort(404, 'File not found')

    @api_route("/files", tags=["client"],
               summary="List the files registered for download",
               responses={200: ("Mapping from registered name to the"
                                " path it resolves to on the server",
                                {"type": "object",
                                 "additionalProperties":
                                     {"type": "string"}})})
    def api_list_all_files(self):
        d = self.get_registered_filenames()

        d = {} if d is None else dict(d)

        return flask.json.jsonify(d), 200

    @api_route("/upload", methods=["POST"], tags=["client"],
               summary="Upload the results of a finished workunit",
               description="Multipart form upload. Besides the files"
                           " themselves, the client sends clientid,"
                           " WUid, a fileinfo JSON object describing"
                           " each file, and, on failure, errorcode and"
                           " failedcommand.",
               responses={200: "Upload accepted",
                          400: "Missing WUid, clientid or fileinfo",
                          403: "A file of that name already exists"})
    def api_upload_file(self):
        if self.read_only:
            flask.abort(409, "This server is only serving the monitoring"
                             " interface; no computation is running")
        clientid = flask.request.form.get('clientid')
        wuid = flask.request.form.get('WUid')
        errorcode = flask.request.form.get('errorcode')
        failedcommand = flask.request.form.get('failedcommand')

        if clientid is None or wuid is None:
            flask.abort(400, "missing WUid and/or clientid")

        # XXX uploaded_files should be a list of tuples.
        # "filename", "path", "type", "command"

        fileinfo = json.loads(flask.request.form.get('fileinfo', "{}"))

        uploaded_files = []
        for fkey, f in flask.request.files.items():
            filename = secure_filename(f.filename)
            path = os.path.join(self.get_upload_folder(filename),
                                filename)
            fi = fileinfo.get(filename)
            if fi is None:
                self.logger.error("Incomplete answer from client:"
                                  " files=%s, data=%s",
                                  flask.request.files.keys(),
                                  json.dumps(fileinfo, indent=4))
                flask.abort(400, f"missing fileinfo for file {filename}")

            if os.path.isfile(path):
                self.logger.error(f"denied request from {clientid}"
                                  f" to overwrite existing file {path}")
                flask.abort(403, "File already exists")
            f.save(path)
            tup = [filename, path, fi["key"]]
            c = fi.get('command')
            if c is not None:
                tup.append(c)
            uploaded_files.append(tuple(tup))
            # self.logger.info(f"saved result file {path}"
            #                  f" for {wuid} (source={clientid})")

        try:
            self.get_wuaccess().result(wuid,
                                       clientid,
                                       uploaded_files,
                                       errorcode,
                                       failedcommand)
        except wudb.StatusUpdateError:
            self.logger.warning(f'Workunit {wuid} was not currently assigned')
        else:
            self.logger.debug(f'Workunit {wuid} completed,'
                              f' uploaded files {uploaded_files}')

        resp = {'message': 'upload completed'}

        return flask.json.jsonify(resp), 200

    @api_route("/WUstatus/<wuid>", tags=["client"],
               summary="Status of one workunit",
               description="Returns the numeric status of a workunit."
                           " The names of the values are, in order:"
                           " AVAILABLE, ASSIGNED, NEED_RESUBMIT,"
                           " RECEIVED_OK, RECEIVED_ERROR, VERIFIED_OK,"
                           " VERIFIED_ERROR, CANCELLED.",
               responses={
                   200: ("The workunit status",
                         {"type": "object",
                          "properties": {
                              "status": {"type": "integer",
                                         "minimum": 0, "maximum": 7}}}),
                   404: "No such workunit"})
    def api_wu_status(self, wuid):
        res = self.get_wuaccess().query(limit=1, eq={"wuid": wuid})

        if not res:
            flask.abort(404, "wuid does not exist")
        else:
            return flask.json.jsonify({"status": res[0]["status"]}), 200

    # ---- the OpenAPI document, and the web ui that consumes the api ----
    #
    # Neither is token-gated. The document describes the interface and
    # carries nothing about the computation; you need to be able to read
    # it in order to write a client at all. The ui files are static
    # assets which must load before the page can so much as ask for a
    # token. Both are still behind the ip filter.

    @api_route(API_PREFIX + "/openapi.json", tags=["docs"],
               summary="OpenAPI 3.1 description of this api",
               responses={200: ("The OpenAPI document",
                                {"type": "object"})})
    def api_openapi(self):
        return flask.json.jsonify(self.openapi_spec()), 200

    @api_route("/api/docs", tags=["docs"],
               summary="Human-readable rendering of the api description",
               responses={200: "An HTML page"})
    def api_docs(self):
        return flask.send_from_directory(UI_DIRECTORY, "docs.html")

    @api_route("/ui/", tags=["ui"],
               summary="Web dashboard",
               responses={200: "An HTML page"})
    def api_ui_index(self):
        return flask.send_from_directory(UI_DIRECTORY, "index.html")

    @api_route("/ui/<path:path>", tags=["ui"],
               summary="Static assets of the web dashboard",
               responses={200: "The asset", 404: "No such asset"})
    def api_ui_asset(self, path):
        return flask.send_from_directory(UI_DIRECTORY, path)

    def openapi_spec(self):
        """
        Assemble the OpenAPI document for the routes we registered.
        """
        doc = spec.build_spec(self.api_routes,
                              servers=[getattr(self, 'url', None)]
                              if getattr(self, 'url', None) else None,
                              **SPEC_INFO)
        doc["components"].setdefault("schemas", {})
        doc["components"]["schemas"].update(spec_schemas())
        return doc
