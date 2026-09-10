# The cado-nfs api

`cado-nfs.py` runs an HTTP server for the whole duration of a
computation. It serves two audiences that have almost nothing in common,
and the split between them is the main thing to understand here.

**Clients.** `cado-nfs-client.py` fetches workunits, downloads the files
they refer to, and uploads results. Those endpoints are unauthenticated,
because a client is handed a URL and a certificate fingerprint and never
a secret. They are gated by `server.whitelist` alone, exactly as they
always were.

**Whoever is running the computation.** Everything under `/api/v1`, and
the web dashboard at `/ui/` that is built on it, reports where the
computation stands and lets you nudge it. Those endpoints require a
bearer token, and are gated by `server.ui_whitelist` as well.

## Reading the api

The full description is generated from the code and served by the
running server:

    <server>/api/v1/openapi.json     OpenAPI 3.1, no token needed
    <server>/api/docs                the same thing, rendered

The document is assembled at run time by `spec.py` from metadata that
each endpoint carries, with no third-party tooling at all — no apispec,
no marshmallow, no flasgger. It therefore answers on any installation
that can run the server. Feed it to whatever OpenAPI tool you prefer.

Routes and documentation come from a single declaration:

```python
@api_route("/api/v1/clients", tags=["monitoring"], auth=True,
           summary="Clients seen by the server",
           responses={200: ("Client list", CLIENT_LIST_SCHEMA)})
@require_token
def api_clients(self):
    ...
```

so an endpoint cannot be reachable and undocumented. The only way to
drift would be to register a route by hand behind that mechanism's back,
and `tests/scripts/cadofactor/test_api.py` fails if anyone does.

## The token

The server writes a token to `<workdir>/<name>.api-token` with mode
`0600` and reuses it across restarts — but only while it is still a
regular file, owned by you, that nobody else can read. If those
conditions stop holding, the server says so and refuses the file rather
than pretending it still protects anything.

Send it as `Authorization: Bearer <token>`.

When it starts, `cado-nfs.py` logs a ready-to-click link:

    Web UI: https://host:8001/ui/#token=b7f3...c19a

The token travels in the URL *fragment*, which browsers never send: it
appears in no request line, no access log and no `Referer` header. The
page moves it straight into `sessionStorage` and clears the address bar.
If you would rather not have anything token-shaped in the cado-nfs log
either, open `/ui/` and paste the token by hand; the page asks for it.

## Reaching the dashboard

`server.ui_whitelist` defaults to `127.0.0.1/32,::1/128`. The intended
way in from elsewhere is therefore an ssh tunnel:

    ssh -L 8001:localhost:8001 the-server-host

and then open <http://localhost:8001/ui/>. This keeps the client
whitelist — which may well be a whole cluster — from also being the list
of people who can administer the computation.

If you would rather expose it directly, widen it explicitly:

    cado-nfs.py ... server.ui_whitelist=127.0.0.1/32,192.168.0.0/24

## What the actions do, and what they deliberately do not

The interesting action is reclaiming the workunits held by a client that
has gone away, so that the work is handed to somebody else instead of
waiting out `tasks.wutimeout`.

None of these actions resubmits anything itself. `ClientServerTask`
keeps its accounting — `wu_submitted`, `wu_received`, `wu_timedout` — in
a *cached* DB-backed dictionary that the running task owns, and asserts
that the number of outstanding workunits never goes negative. Anything
here that wrote those rows would be invisible to the task, would then be
overwritten by it, and would eventually trip that assertion.

So the actions set exactly one thing: the workunit's status becomes
`NEED_RESUBMIT`. `ClientServerTask.resubmit_timed_out_wus()` already
looks for that on every `wait()` tick and puts the workunit back through
the task's own `cancel_wu()` and `resubmit_one_wu()`, which keep the
counters straight. The effect is therefore not instantaneous: it lands
within `tasks.wutimeoutcheck` seconds, 60 by default, and the api says
so in its answer.

Two things are missing on purpose:

* **There is no raw cancel.** Setting `CANCELLED` from outside would
  never decrement the task's `wu_timedout`, so
  `get_number_outstanding_wus()` would never reach zero and the task
  would wait forever.

* **Resubmitting a workunit that belongs to a task other than the one
  running is refused with 409.** `resubmit_timed_out_wus()` queries
  `NEED_RESUBMIT` across the whole table without filtering by task, and
  charges whatever it finds to the task that is running. Marking a
  leftover from an earlier phase would corrupt that task's accounting.
  Such leftovers do occur: a `SievingTask` that finishes early because
  it already has enough relations drops its outstanding workunits
  without cancelling them.

## Where things live

    api/spec.py     the @api_route decorator, and OpenAPI assembly
    api/auth.py     the token, and the require_token decorator
    api/views.py    read-only aggregation over the database
    api/admin.py    the /api/v1 endpoints
    api/ui/         the web dashboard (no build step, no vendored code)
    ../api_server.py  the flask application, TLS and WSGI plumbing

`views.py` is where the database quirks are handled, and its module
docstring lists them: attaching to a state dictionary with
`DictDbDirectAccess` takes an EXCLUSIVE lock, the MySQL cursor rewrites
`ASC` and `purge` on their way to the server, and the parameter
placeholder differs between backends. Read it before adding a query.

## Dependencies

| Component | Licence | Needed for |
|---|---|---|
| flask, werkzeug, jinja2 | BSD-3-Clause | the server at all |
| requests | Apache-2.0 | `cado-nfs-client.py` |
| rich | MIT | the live view of `cado-nfs-monitor.py` only |

Nothing else. The dashboard is plain ES modules and CSS that are part of
cado-nfs; there is no package manager, no build step and no vendored
third-party JavaScript, and the charts are inline SVG drawn by hand.
`cado-nfs-monitor.py` needs nothing beyond the standard library and can
be copied to another machine on its own.
