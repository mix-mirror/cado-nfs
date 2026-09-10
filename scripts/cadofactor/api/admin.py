"""
The authenticated part of the cado-nfs api, under /api/v1.

These endpoints are what the web ui at /ui/ and cado-nfs-monitor.py are
built on. There is deliberately no second path: everything either of
them displays comes through here, so that the api cannot quietly rot
while the ui keeps working.

A word on the actions, because they are the delicate part.

ClientServerTask keeps its accounting -- wu_submitted, wu_received,
wu_timedout -- in a *cached* DB-backed dictionary, and asserts that the
number of outstanding workunits never goes negative. Anything here that
wrote those rows would be invisible to the running task, would then be
overwritten by it, and would eventually trip that assertion. So this
module never touches task state, and never resubmits anything itself.

It does not need to. cadotask.ClientServerTask.resubmit_timed_out_wus()
already looks for workunits in state NEED_RESUBMIT on every wait() tick,
and puts them back through the task's own cancel_wu() and
resubmit_one_wu(), which do keep the counters straight. Setting that one
status is the whole of what the actions below do.

The same method queries NEED_RESUBMIT across the *entire* table without
filtering by task, and then charges the work to whichever task is
running. Marking a workunit left over from an earlier phase would
therefore corrupt the current task's accounting, so every action here
refuses workunits that do not belong to the task currently running.
"""

import hashlib
import json
import logging
import os
import time

import flask

from cadofactor.api.auth import require_token
from cadofactor.api.spec import api_route, json_body, query_parameter
from cadofactor.api.views import (TOP_CLIENTS_LIMIT,
                                  TUNABLE_PARAMETERS)
from cadofactor.workunit import STATUS_NAMES, WuStatus

logger = logging.getLogger("API server")

API = "/api/v1"

# Bound on how much of the log file one request may pull.
MAX_LOG_LINES = 5000

# Bound on the captured output shown with a single workunit. This is
# the tail of what a client's command printed, which is the first thing
# anybody wants when a workunit failed.
WU_OUTPUT_LINES = 60

# Bound on a single page of workunits.
MAX_PAGE = 500

# Bound on a single page of clients, and how many of them the summary
# form names. Both exist because pools get big: 1400 clients is half a
# megabyte of JSON, and that is a small run.
MAX_CLIENT_PAGE = 500
TOP_CLIENTS = TOP_CLIENTS_LIMIT


def spec_schemas():
    """
    Component schemas referenced from several operations.
    """
    return {
        "Workunit": {
            "type": "object",
            "description": "A unit of work as handed to a client. The"
                           " canonical description lives in"
                           " scripts/cadofactor/workunit.py.",
            "required": ["id"],
            "properties": {
                "id": {"type": "string"},
                "commands": {"type": "array",
                             "items": {"type": "string"}},
                "deadline": {"type": "number"},
                "files": {"type": "object"},
            },
        },
        "Error": {
            "type": "object",
            "properties": {
                "code": {"type": "integer"},
                "name": {"type": "string"},
                "description": {"type": "string"},
            },
        },
        "TaskProgress": {
            "type": "object",
            "properties": {
                "name": {"type": "string"},
                "title": {"type": "string"},
                "phase": {"type": "string",
                          "enum": ["done", "running", "pending"]},
                "achievement": {"type": "number",
                                "description": "Fraction in [0,1], as"
                                               " the task itself"
                                               " computes it"},
                "eta": {"type": "string"},
                "highlights": {"type": "object"},
                "times": {"type": "object"},
                "stats": {"type": "array", "items": {"type": "string"}},
            },
        },
        "Client": {
            "type": "object",
            "properties": {
                "clientid": {"type": "string"},
                "in_flight": {"type": "integer"},
                "completed": {"type": "integer"},
                "failed": {"type": "integer"},
                "share": {"type": "number"},
                "last_seen": {"type": ["number", "null"],
                              "description": "epoch seconds; subtract"
                                             " from the"
                                             " X-Cado-Server-Time"
                                             " response header to get"
                                             " an age"},
                "oldest_assignment": {"type": ["number", "null"],
                                      "description": "epoch seconds of"
                                                     " the longest"
                                                     " outstanding"
                                                     " assignment"},
                "state": {"type": "string",
                          "enum": ["working", "idle", "stale", "gone",
                                   "unknown"]},
                "typical_turnaround": {
                    "type": ["number", "null"],
                    "description": "Median seconds this client has"
                                   " recently taken to return a"
                                   " workunit, or null if too few"
                                   " samples"},
                "turnaround_samples": {
                    "type": "integer",
                    "description": "How many recent workunits that"
                                   " median rests on"},
                "stale_after": {
                    "type": "number",
                    "description": "Seconds of silence after which this"
                                   " client is called stale. Derived"
                                   " from its own turnaround when known,"
                                   " capped by tasks.wutimeout."},
                "liveness_basis": {
                    "type": "string",
                    "enum": ["turnaround", "wutimeout"],
                    "description": "Which of the two the threshold came"
                                   " from"},
            },
        },
        "WorkunitInfo": {
            "type": "object",
            "properties": {
                "wuid": {"type": "string"},
                "status": {"type": "integer"},
                "status_name": {"type": "string"},
                "task": {"type": ["string", "null"]},
                "identifier": {"type": ["string", "null"]},
                "attempt": {"type": "integer"},
                "assignedclient": {"type": ["string", "null"]},
                "resultclient": {"type": ["string", "null"]},
                "timecreated": {"type": ["number", "null"]},
                "timeassigned": {"type": ["number", "null"]},
                "timeresult": {"type": ["number", "null"]},
            },
        },
        "ActionResult": {
            "type": "object",
            "properties": {
                "marked": {"type": "array", "items": {"type": "string"}},
                "skipped": {"type": "array",
                            "items": {"type": "object"}},
                "message": {"type": "string"},
            },
        },
    }


def next_snapshot(workdir, name, overrides):
    """
    Record the parameters now in force as a fresh snapshot.

    cado-nfs runs are meant to be reproducible from a
    <name>.parameters_snapshot.<N> file, and cado-nfs.py writes one at
    the start of every run. A ceiling raised from the api would
    otherwise appear in no snapshot at all, so we write the next one in
    the sequence: the highest-numbered snapshot keeps describing the
    parameters actually in force, and resuming from it reproduces this
    configuration rather than the one the run started with.

    Returns the path written, or None if there was no snapshot to build
    on -- in which case the caller should say so rather than pretend
    the record is complete.
    """
    from cadofactor import cadoparams

    existing = []
    prefix = name + ".parameters_snapshot."
    try:
        for entry in os.listdir(workdir):
            if entry.startswith(prefix) and entry[len(prefix):].isdigit():
                existing.append(int(entry[len(prefix):]))
    except OSError as e:
        logger.error("Cannot list %s to write a snapshot: %s", workdir, e)
        return None
    if not existing:
        return None

    latest = os.path.join(workdir, prefix + str(max(existing)))
    parameters = cadoparams.Parameters()
    try:
        parameters.readfile(latest)
    except Exception as e:
        logger.error("Cannot read %s: %s", latest, e)
        return None

    for key, value in overrides.items():
        parameters.set_simple("tasks." + key, value)

    # O_EXCL, and step forward on a clash: cado-nfs.py picks the first
    # free number the same way, and it may be doing so right now.
    for number in range(max(existing) + 1, max(existing) + 100):
        path = os.path.join(workdir, prefix + str(number))
        try:
            fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o666)
        except FileExistsError:
            continue
        except OSError as e:
            logger.error("Cannot write %s: %s", path, e)
            return None
        with os.fdopen(fd, "w") as handle:
            handle.write(str(parameters))
            handle.write("\n")
        return path
    return None


def tail_file(path, lines):
    """
    Return the last `lines` lines of a file, without reading it whole.

    Log files of a long computation reach hundreds of megabytes, so we
    walk backwards in blocks rather than doing readlines().

    >>> import tempfile, os
    >>> fd, p = tempfile.mkstemp()
    >>> _ = os.write(fd, b''.join(b'line %d\\n' % i for i in range(1000)))
    >>> os.close(fd)
    >>> t = tail_file(p, 3)
    >>> t
    ['line 997', 'line 998', 'line 999']
    >>> len(tail_file(p, 10000))
    1000
    >>> tail_file(p, 0)
    []
    >>> tail_file('/nonexistent/at/all', 5)
    []
    >>> os.unlink(p)
    """
    if lines <= 0:
        return []
    block = 8192
    try:
        with open(path, "rb") as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            data = b""
            while size > 0 and data.count(b"\n") <= lines:
                step = min(block, size)
                size -= step
                f.seek(size)
                data = f.read(step) + data
    except OSError:
        return []
    text = data.decode("utf-8", "replace")
    return text.splitlines()[-lines:]


def json_response(payload, status=200):
    """
    Serialize a payload, with an ETag so that polling is cheap.

    The dashboard refreshes a few panels every couple of seconds; most
    of those refreshes find nothing changed, and should cost a 304
    rather than a re-render.

    For that to work at all, response bodies must not contain anything
    that changes on every single call. So they carry absolute
    timestamps, never elapsed times, and the server's own clock travels
    in a header which is deliberately left out of the ETag. Consumers
    subtract the two, which also spares them any clock skew between the
    server and whoever is looking at the dashboard.
    """
    body = json.dumps(payload, sort_keys=True, default=str)
    etag = hashlib.sha1(body.encode("utf-8")).hexdigest()
    # Go through werkzeug's own accessors rather than comparing the raw
    # header: an ETag travels quoted on the wire, and set_etag /
    # if_none_match are what get the quoting right in both directions.
    if flask.request.if_none_match.contains(etag):
        response = flask.make_response("", 304)
    else:
        response = flask.make_response(body, status)
        response.content_type = "application/json"
    response.set_etag(etag)
    response.headers["Cache-Control"] = "no-cache"
    # Set after set_etag, and never folded into the body, so that an
    # otherwise unchanged answer keeps revalidating as a 304.
    response.headers["X-Cado-Server-Time"] = repr(time.time())
    return response


def _int_arg(name, default, minimum=None, maximum=None):
    """
    Read a bounded integer from the query string.
    """
    raw = flask.request.args.get(name)
    if raw is None or raw == "":
        return default
    try:
        value = int(raw)
    except ValueError:
        flask.abort(400, "%s must be an integer" % name)
    if minimum is not None and value < minimum:
        value = minimum
    if maximum is not None and value > maximum:
        value = maximum
    return value


def _status_arg():
    """
    Read a status filter given either by name or by number.
    """
    raw = flask.request.args.get("status")
    if raw is None or raw == "":
        return None
    if raw.isdigit():
        value = int(raw)
    else:
        try:
            value = STATUS_NAMES.index(raw.upper())
        except ValueError:
            flask.abort(400, "unknown status %r; expected one of %s"
                        % (raw, ", ".join(STATUS_NAMES)))
    if not 0 <= value < len(STATUS_NAMES):
        flask.abort(400, "status out of range")
    return value


class AdminEndpoints(object):
    """
    Monitoring and administration endpoints, bound to an ApiServer.

    They live on a separate object purely for readability; they are
    registered on the same flask application, through the same
    @api_route mechanism as the client endpoints.
    """

    def __init__(self, app):
        self.app = app

    # ---------------- helpers ----------------

    @property
    def views(self):
        return self.app.views

    def _current_task(self):
        return self.views.progress_state().get("current") or None

    def _captured_output(self, files):
        """
        The tail of whatever the client's commands printed.

        Paths come out of the database, where the server itself put
        them when it saved the upload, so they should already be inside
        the working directory -- but this reads files off disk on
        request, so it checks rather than assumes.
        """
        root = os.path.realpath(self.app.wdir) if self.app.wdir else None
        out = []
        for entry in files or []:
            kind = (entry.get("type") or "")
            if not (kind.startswith("stderr") or kind.startswith("stdout")):
                continue
            path = entry.get("path")
            if not path:
                continue
            if root is not None:
                real = os.path.realpath(path)
                if not (real == root or real.startswith(root + os.sep)):
                    logger.warning("refusing to read %s: outside %s",
                                   path, root)
                    continue
            out.append({"type": kind,
                        "filename": entry.get("filename"),
                        "lines": tail_file(path, WU_OUTPUT_LINES)})
        return out

    def _refuse_if_read_only(self):
        """
        Actions need a running task to pick them up.

        In --ui-only mode there is none, so a NEED_RESUBMIT we set here
        would simply sit there until somebody started a computation
        again -- at which point it would be charged to whatever task
        happened to be running. Refuse instead of leaving that trap.
        """
        if getattr(self.app, "read_only", False):
            flask.abort(409, "This server is only serving the monitoring"
                             " interface; no task is running to act on"
                             " what you would be asking for")

    def _mark_for_resubmit(self, rows):
        """
        Put workunits back in the pool, by way of NEED_RESUBMIT.

        Returns (marked, skipped). A workunit is skipped when it is not
        currently assigned, or when it belongs to a task other than the
        one running -- see this module's docstring for why the latter
        matters.
        """
        current = self._current_task()
        wuaccess = self.app.get_wuaccess()
        marked, skipped = [], []
        for row in rows:
            wuid = row["wuid"]
            if row["status"] != WuStatus.ASSIGNED:
                skipped.append({"wuid": wuid,
                                "reason": "not assigned (status is %s)"
                                          % row["status_name"]})
                continue
            if current is None:
                skipped.append({"wuid": wuid,
                                "reason": "no task is running, so nothing"
                                          " would pick this up"})
                continue
            if row["task"] != current:
                skipped.append({"wuid": wuid,
                                "reason": "belongs to task %s, but %s is"
                                          " running; resubmitting it would"
                                          " be charged to the wrong task"
                                          % (row["task"], current)})
                continue
            wuaccess.set_status(WuStatus.NEED_RESUBMIT,
                                eq={"wuid": wuid,
                                    "status": WuStatus.ASSIGNED})
            marked.append(wuid)
        if marked:
            self.views.invalidate()
            logger.info("api: marked %d workunit(s) for resubmission: %s",
                        len(marked), ", ".join(marked))
        return marked, skipped

    def _action_result(self, marked, skipped):
        if marked:
            message = ("%d workunit(s) marked for resubmission; the"
                       " running task picks them up on its next timeout"
                       " check." % len(marked))
        else:
            message = "Nothing was marked for resubmission."
        return json_response({"marked": marked,
                              "skipped": skipped,
                              "message": message})

    # ---------------- monitoring ----------------

    @api_route(API + "/info", tags=["monitoring"], auth=True,
               summary="What computation this server belongs to",
               responses={200: ("Identity and top-level timings",
                                {"type": "object"})})
    @require_token
    def api_info(self):
        state = self.views.read_state_table("tasks")
        published = self.views.progress_state()
        started = state.get("starttime")
        return json_response({
            "name": self.app.computation_name,
            "workdir": self.app.wdir,
            "url": getattr(self.app, "url", None),
            "database": self.app.database_uri.uri_without_credentials,
            "computation": published.get("computation"),
            "computation_desc": published.get("computation_desc"),
            "algorithm": published.get("algo"),
            "N": state.get("N"),
            "starttime": started,
            "elapsed_before": state.get("elapsed"),
            "serving_workunits": self.app.serving.get(),
            "wutimeout": self.views.wutimeout(),
            "api_version": 1,
        })

    @api_route(API + "/progress", tags=["monitoring"], auth=True,
               summary="Where the computation stands",
               description="The task pipeline in order, which task is"
                           " running, and how far along it is. The"
                           " achievement and eta fields are the very"
                           " numbers the task itself logs.",
               responses={200: ("Pipeline progress",
                                {"type": "object",
                                 "properties": {
                                     "current": {"type": ["string",
                                                          "null"]},
                                     "finished": {"type": "boolean"},
                                     "tasks": {
                                         "type": "array",
                                         "items": {"$ref": "#/components"
                                                           "/schemas"
                                                           "/TaskProgress"}
                                     }}})})
    @require_token
    def api_progress(self):
        return json_response(self.views.progress())

    @api_route(API + "/workunits/summary", tags=["monitoring"], auth=True,
               summary="How many workunits are in each status",
               responses={200: ("Counts keyed by status name",
                                {"type": "object",
                                 "additionalProperties":
                                     {"type": "integer"}})})
    @require_token
    def api_workunits_summary(self):
        counts = self.views.status_counts()
        return json_response({
            "counts": counts,
            "total": sum(counts.values()),
            "outstanding": (counts.get("AVAILABLE", 0)
                            + counts.get("ASSIGNED", 0)
                            + counts.get("NEED_RESUBMIT", 0)),
        })

    @api_route(API + "/workunits", tags=["monitoring"], auth=True,
               summary="Browse the workunits table",
               parameters=[
                   query_parameter("status", {"type": "string"},
                                   "Status name or number to filter on"),
                   query_parameter("assigned_to", {"type": "string"},
                                   "Only workunits handed to this"
                                   " client and not yet returned"),
                   query_parameter("result_from", {"type": "string"},
                                   "Only workunits returned by this"
                                   " client"),
                   query_parameter("task", {"type": "string"},
                                   "Only workunits of this task"),
                   query_parameter("assigned_older_than",
                                   {"type": "integer"},
                                   "Only workunits assigned more than"
                                   " this many seconds ago"),
                   query_parameter("limit", {"type": "integer"},
                                   "Page size, at most %d" % MAX_PAGE),
                   query_parameter("offset", {"type": "integer"},
                                   "Rows to skip"),
               ],
               responses={200: ("A page of workunits",
                                {"type": "object",
                                 "properties": {
                                     "workunits": {
                                         "type": "array",
                                         "items": {"$ref": "#/components"
                                                           "/schemas"
                                                           "/WorkunitInfo"}
                                     }}}),
                          400: "Bad filter value"})
    @require_token
    def api_workunits(self):
        limit = _int_arg("limit", 50, minimum=1, maximum=MAX_PAGE)
        offset = _int_arg("offset", 0, minimum=0)
        older = _int_arg("assigned_older_than", None, minimum=0)
        rows = self.views.list_workunits(
            status=_status_arg(),
            assigned_to=flask.request.args.get("assigned_to") or None,
            result_from=flask.request.args.get("result_from") or None,
            task=flask.request.args.get("task") or None,
            assigned_older_than=older,
            limit=limit,
            offset=offset)
        return json_response({"workunits": rows,
                              "limit": limit,
                              "offset": offset,
                              "returned": len(rows)})

    @api_route(API + "/workunits/<wuid>", tags=["monitoring"], auth=True,
               summary="One workunit in full",
               description="The workunit body, its files, every attempt"
                           " that has been made at the same piece of"
                           " work, and the tail of whatever the"
                           " client's commands printed -- which is"
                           " normally the first thing wanted when one"
                           " has failed.",
               parameters=[
                   query_parameter("output", {"type": "boolean"},
                                   "Include the captured stdout and"
                                   " stderr (default: yes)"),
               ],
               responses={
                   200: ("The workunit",
                         {"$ref": "#/components/schemas/WorkunitInfo"}),
                   404: "No such workunit"})
    @require_token
    def api_workunit(self, wuid):
        rows = self.app.get_wuaccess().query(limit=1, eq={"wuid": wuid})
        if not rows:
            flask.abort(404, "wuid does not exist")
        detail = self.views.workunit_detail(rows[0])
        if flask.request.args.get("output", "1") not in ("0", "false",
                                                         "no"):
            detail["output"] = self._captured_output(detail.get("files"))
        return json_response(detail)

    @api_route(API + "/clients/<clientid>", tags=["monitoring"],
               auth=True,
               summary="One client in full",
               description="Everything /api/v1/clients reports about"
                           " this one, plus the workunits it is holding"
                           " and a page of what it has recently handed"
                           " back, with how long each took. Bounded by"
                           " construction, which is why it can afford"
                           " to be detailed where the list cannot.",
               responses={
                   200: ("The client",
                         {"$ref": "#/components/schemas/Client"}),
                   404: "No such client has been seen"})
    @require_token
    def api_client(self, clientid):
        detail = self.views.client_detail(clientid)
        if detail is None:
            flask.abort(404, "no client called %r has asked for work"
                        % clientid)
        return json_response(detail)

    @api_route(API + "/clients", tags=["monitoring"], auth=True,
               summary="Clients seen by the server",
               description="Derived from the workunits table: what each"
                           " client holds now, what it has returned,"
                           " when we last heard from it, and whether it"
                           " still looks alive."
                           " The staleness threshold is per client:"
                           " a client that has returned enough"
                           " workunits has shown how long its workunits"
                           " take, and being silent for several times"
                           " that is suspicious well before"
                           " tasks.wutimeout would say so. See"
                           " stale_after and liveness_basis. The"
                           " threshold is never longer than"
                           " tasks.wutimeout, since the work is"
                           " reassigned at that point regardless, and"
                           " while few workunits back the estimate it"
                           " is also floored at a small multiple of"
                           " tasks.wutimeoutcheck."
                           " Pools get large -- 1400 clients is half a"
                           " megabyte of JSON -- so the full list is"
                           " paginated, and summary=1 returns only the"
                           " tallies and the top contributors, which is"
                           " what a dashboard should poll.",
               parameters=[
                   query_parameter("summary", {"type": "boolean"},
                                   "Return only counts by state and the"
                                   " top %d contributors, omitting the"
                                   " per-client array" % TOP_CLIENTS),
                   query_parameter("limit", {"type": "integer"},
                                   "Page size, at most %d"
                                   % MAX_CLIENT_PAGE),
                   query_parameter("offset", {"type": "integer"},
                                   "Clients to skip"),
               ],
               responses={200: ("Client list",
                                {"type": "object",
                                 "properties": {
                                     "clients": {
                                         "type": "array",
                                         "items": {"$ref": "#/components"
                                                           "/schemas"
                                                           "/Client"}
                                     }}})})
    @require_token
    def api_clients(self):
        clients = self.views.clients()
        by_state = {}
        for client in clients:
            by_state[client["state"]] = by_state.get(client["state"], 0) + 1
        payload = {
            "counts": by_state,
            "total": len(clients),
            "wutimeout": self.views.wutimeout(),
        }

        # A real pool is large: a c180 polyselect run with 1400 clients
        # makes this endpoint half a megabyte, and a dashboard that
        # asked for it every couple of seconds would be taking capacity
        # away from the computation it is supposed to be watching. So
        # the whole list is opt-in, and paginated.
        if flask.request.args.get("summary") in ("1", "true", "yes"):
            # Ten rows, so a couple of extra fields cost nothing and
            # save the caller from needing the full list to draw a
            # useful table.
            payload["top"] = [
                {k: c[k] for k in ("clientid", "state", "completed",
                                   "failed", "in_flight", "share",
                                   "last_seen", "typical_turnaround")}
                for c in clients[:TOP_CLIENTS]]
            payload["truncated"] = max(0, len(clients) - TOP_CLIENTS)
            return json_response(payload)

        limit = _int_arg("limit", MAX_CLIENT_PAGE, minimum=1,
                         maximum=MAX_CLIENT_PAGE)
        offset = _int_arg("offset", 0, minimum=0)
        payload["clients"] = clients[offset:offset + limit]
        payload["limit"] = limit
        payload["offset"] = offset
        payload["returned"] = len(payload["clients"])
        return json_response(payload)

    @api_route(API + "/stats", tags=["monitoring"], auth=True,
               summary="Per-task cpu and elapsed times, and the"
                       " statistics each task reports about itself",
               responses={200: ("Statistics per task",
                                {"type": "object"})})
    @require_token
    def api_stats(self):
        progress = self.views.progress()
        out = []
        for task in progress["tasks"]:
            entry = {"name": task["name"], "title": task["title"]}
            if "times" in task:
                entry["times"] = task["times"]
            if "stats" in task:
                entry["stats"] = task["stats"]
            out.append(entry)
        return json_response({"tasks": out,
                              "elapsed": progress.get("elapsed"),
                              "cputotal": progress.get("cputotal")})

    @api_route(API + "/log", tags=["monitoring"], auth=True,
               summary="Tail of the computation's log file",
               parameters=[
                   query_parameter("tail", {"type": "integer"},
                                   "Number of trailing lines, at most"
                                   " %d" % MAX_LOG_LINES),
               ],
               responses={200: ("Log lines",
                                {"type": "object",
                                 "properties": {
                                     "lines": {"type": "array",
                                               "items": {
                                                   "type": "string"}}}}),
                          404: "The server does not know where its log is"})
    @require_token
    def api_log(self):
        if not self.app.wdir or not self.app.computation_name:
            flask.abort(404, "no log file is associated with this server")
        path = os.path.join(self.app.wdir,
                            self.app.computation_name + ".log")
        count = _int_arg("tail", 200, minimum=1, maximum=MAX_LOG_LINES)
        return json_response({"path": path,
                              "lines": tail_file(path, count)})

    # ---------------- actions ----------------

    @api_route(API + "/workunits/<wuid>/resubmit", methods=["POST"],
               tags=["actions"], auth=True,
               summary="Put one workunit back in the pool",
               description="Marks the workunit NEED_RESUBMIT. The"
                           " running task notices on its next timeout"
                           " check -- within tasks.wutimeoutcheck"
                           " seconds, 60 by default -- and resubmits it"
                           " itself, so that its accounting stays"
                           " correct.",
               responses={
                   200: ("What was marked",
                         {"$ref": "#/components/schemas/ActionResult"}),
                   404: "No such workunit",
                   409: "The workunit is not assigned, or belongs to a"
                        " task that is not currently running"})
    @require_token
    def api_workunit_resubmit(self, wuid):
        self._refuse_if_read_only()
        rows = self.app.get_wuaccess().query(limit=1, eq={"wuid": wuid})
        if not rows:
            flask.abort(404, "wuid does not exist")
        described = [self.views.describe_workunit(rows[0])]
        marked, skipped = self._mark_for_resubmit(described)
        if not marked:
            flask.abort(409, skipped[0]["reason"])
        return self._action_result(marked, skipped)

    @api_route(API + "/workunits/reclaim", methods=["POST"],
               tags=["actions"], auth=True,
               summary="Reclaim workunits that have been out too long",
               description="Marks every workunit assigned longer ago"
                           " than older_than seconds NEED_RESUBMIT. This"
                           " is the bulk form of the per-client reclaim,"
                           " for when several machines went away at"
                           " once.",
               request_body=json_body(
                   {"type": "object",
                    "properties": {
                        "older_than": {
                            "type": "integer",
                            "description": "Age in seconds of the"
                                           " assignment. Defaults to"
                                           " tasks.wutimeout."}}},
                   required=False),
               responses={200: ("What was marked",
                                {"$ref": "#/components/schemas"
                                         "/ActionResult"})})
    @require_token
    def api_workunits_reclaim(self):
        self._refuse_if_read_only()
        body = flask.request.get_json(silent=True) or {}
        try:
            older = float(body.get("older_than")
                          or self.views.wutimeout())
        except (TypeError, ValueError):
            flask.abort(400, "older_than must be a number of seconds")
        rows = self.views.list_workunits(status=WuStatus.ASSIGNED,
                                         assigned_older_than=older,
                                         limit=MAX_PAGE)
        marked, skipped = self._mark_for_resubmit(rows)
        return self._action_result(marked, skipped)

    @api_route(API + "/clients/<clientid>/reclaim", methods=["POST"],
               tags=["actions"], auth=True,
               summary="Reclaim the workunits a client is holding",
               description="For a client that has gone away. Marks"
                           " everything currently assigned to it"
                           " NEED_RESUBMIT, so that the running task"
                           " hands the work to somebody else instead of"
                           " waiting out the full timeout."
                           " If the client turns out to be alive after"
                           " all, its upload of the reclaimed workunit"
                           " is refused with 403 because the"
                           " replacement got there first; it gives up"
                           " on that one and carries on. Nothing is"
                           " lost.",
               responses={200: ("What was marked",
                                {"$ref": "#/components/schemas"
                                         "/ActionResult"})})
    @require_token
    def api_client_reclaim(self, clientid):
        self._refuse_if_read_only()
        rows = self.views.list_workunits(status=WuStatus.ASSIGNED,
                                         assigned_to=clientid,
                                         limit=MAX_PAGE)
        marked, skipped = self._mark_for_resubmit(rows)
        return self._action_result(marked, skipped)

    @api_route(API + "/parameters", tags=["actions"], auth=True,
               summary="The ceilings that may be raised while running",
               description="A deliberately tiny list. Both entries do"
                           " nothing except abort the computation when"
                           " a counter passes them, so raising one"
                           " cannot change what the run computes, only"
                           " whether it survives long enough to finish."
                           " Everything else stays in the parameter"
                           " file, so that a run remains reproducible"
                           " from its snapshot.",
               responses={200: ("Current value, default, and how close"
                                " the running task is to it",
                                {"type": "object"})})
    @require_token
    def api_parameters(self):
        return json_response({"parameters": self.views.tunables()})

    @api_route(API + "/parameters/<name>", methods=["POST"],
               tags=["actions"], auth=True,
               summary="Raise one ceiling",
               description="Writes the new value where the running"
                           " tasks will pick it up, and records it by"
                           " writing the next parameters_snapshot in"
                           " the sequence, so that the highest-numbered"
                           " snapshot still describes the parameters in"
                           " force and resuming from it reproduces"
                           " them. The change is also logged.",
               request_body=json_body(
                   {"type": "object",
                    "required": ["value"],
                    "properties": {
                        "value": {"type": "integer", "minimum": 1}}}),
               responses={
                   200: ("The new value, and the snapshot written",
                         {"type": "object"}),
                   400: "Not an integer, or not a positive one",
                   404: "No such tunable parameter",
                   409: "The value would abort the run immediately, or"
                        " no task is running to pick it up"})
    @require_token
    def api_set_parameter(self, name):
        self._refuse_if_read_only()
        if name not in TUNABLE_PARAMETERS:
            flask.abort(404, "%r is not tunable while the computation"
                             " runs. Tunable: %s. Everything else"
                             " belongs in the parameter file, so that"
                             " the run stays reproducible from its"
                             " snapshot."
                        % (name, ", ".join(sorted(TUNABLE_PARAMETERS))))

        body = flask.request.get_json(silent=True) or {}
        if "value" not in body:
            flask.abort(400, "expected a JSON body with a 'value'")
        try:
            value = int(body["value"])
        except (TypeError, ValueError):
            flask.abort(400, "'value' must be an integer")
        if value < 1:
            flask.abort(400, "'value' must be positive")

        tunables = self.views.tunables()
        counter = tunables[name]["counter_value"]
        if counter is not None and value <= counter:
            # Lowering a ceiling below where the counter already stands
            # would abort the computation the next time it moves. That
            # is almost certainly not what was meant.
            flask.abort(409,
                        "%s already stands at %d, so setting %s to %d"
                        " would abort the computation at the next"
                        " occurrence. Pass a larger value."
                        % (tunables[name]["counter"], counter, name,
                           value))

        overrides = self.views.overrides()
        overrides[name] = value

        # Write the record before the value takes effect, so that a
        # crash in between leaves a snapshot promising something the
        # run did not do, rather than a run doing something no snapshot
        # records.
        snapshot = None
        if self.app.wdir and self.app.computation_name:
            snapshot = next_snapshot(self.app.wdir,
                                     self.app.computation_name,
                                     overrides)
            if snapshot is None:
                logger.warning("Could not write a parameters snapshot"
                               " for the change to %s", name)

        self.app.session().overrides[name] = value
        self.views.invalidate()
        logger.info("api: %s raised to %d (was %d)%s",
                    name, value, tunables[name]["value"],
                    "; recorded in " + os.path.basename(snapshot)
                    if snapshot else "; NO snapshot could be written")

        return json_response({
            "name": name,
            "value": value,
            "previous": tunables[name]["value"],
            "snapshot": snapshot,
            "message": ("%s is now %d. The running task picks this up"
                        " the next time it checks." % (name, value))
            + ("" if snapshot else
               " WARNING: no parameters snapshot could be written, so"
               " this change is recorded only in the log."),
        })

    @api_route(API + "/serving", methods=["GET", "POST"],
               tags=["actions"], auth=True,
               summary="Whether the server hands workunits out",
               description="POST {\"serving\": false} to stop giving"
                           " clients new work; they will be told 410"
                           " and will terminate. This is what cado-nfs"
                           " does itself at the end of the distributed"
                           " phase.",
               request_body=json_body(
                   {"type": "object",
                    "required": ["serving"],
                    "properties": {"serving": {"type": "boolean"}}},
                   required=False),
               responses={200: ("The resulting state",
                                {"type": "object",
                                 "properties": {
                                     "serving": {"type": "boolean"}}}),
                          400: "Missing or non-boolean serving field"})
    @require_token
    def api_serving(self):
        if flask.request.method == "POST":
            self._refuse_if_read_only()
            body = flask.request.get_json(silent=True) or {}
            if "serving" not in body:
                flask.abort(400, "expected a JSON body with a"
                                 " 'serving' boolean")
            wanted = body["serving"]
            if not isinstance(wanted, bool):
                flask.abort(400, "'serving' must be a boolean")
            if wanted:
                self.app.resume_serving_wus()
            else:
                self.app.stop_serving_wus()
        return json_response({"serving": self.app.serving.get()})


if __name__ == "__main__":
    import doctest
    doctest.testmod()
