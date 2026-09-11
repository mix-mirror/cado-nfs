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

### The dashboard's urls

Everything the dashboard is showing is in the fragment, so every view
and every filter is a link that can be bookmarked or sent to somebody:

    #overview
    #stage/sieving                       one phase of the pipeline
    #clients?group_by=cluster            the roll-up
    #group/alpha?group_by=cluster        one cluster, and its mass action
    #clients?state=stale,gone            only the quiet ones
    #clients/alpha-15+101                one client
    #workunits?task=sieving&status=ASSIGNED&client=alpha-15
    #workunits/c60_sieving_100-200       one workunit
    #log?tail=500

The fragment is also where the token arrives (`/ui/#token=...`), which
`app.js` moves into session storage and wipes from the address bar
before any of the above is parsed.

## Deciding that a client has gone away

`tasks.wutimeout` is one number for the whole computation, chosen so
that the slowest machine in the pool is not cheated. As a staleness
signal it is therefore blunt: it says nothing about any particular
client, and on a heterogeneous pool it is far too patient for most of
them.

But a client that has been running for a while has already told us how
long its workunits take — the workunits table records `timeassigned`
and `timeresult` for everything it has handed back. So the threshold is
per client: `stale_after` is `TURNAROUND_FACTOR` (6) times the median
turnaround of that client's recent workunits, capped at
`tasks.wutimeout`, past which the server reassigns the work anyway and
calling the client "working" would be a lie.

While that estimate still rests on fewer than
`TURNAROUND_TRUSTED_SAMPLES` (10) workunits it is also floored, at
`MIN_STALE_AFTER_CHECKS` (2) times `tasks.wutimeoutcheck` — the
interval at which the running task actually looks for work to reassign,
so there is nothing to be gained by being twitchier than a small
multiple of it. Once enough workunits agree with each other the floor
is dropped: consistent evidence is precisely the case where an
arbitrary floor has no rationale left.

A machine that returns a workunit every two minutes and has been silent
for twenty is thus flagged long before a three-hour `wutimeout` would
notice, which is the whole point.

The sample is bounded and recent — the last `TURNAROUND_WINDOW` (2000)
finished workunits overall, walked backwards from the newest by primary
key. Recent matters as much as bounded: turnaround depends on which
task is running, and polyselect and sieving workunits are nothing
alike. The cost is that a client which finishes rarely may have nothing
inside the window; it then falls back to `wutimeout`, and says so.
`/api/v1/clients` reports `typical_turnaround`, `turnaround_samples`,
`stale_after` and `liveness_basis` for every client, so the verdict can
be checked rather than taken on faith — and both uis show the
reasoning.

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

One consequence is worth knowing about. If you reclaim a client that
turns out not to have gone away after all, it will eventually finish
the workunit and try to upload a result file that the replacement has
already produced. The server refuses the overwrite with 403, the client
retries with exponential backoff and then gives up on that upload and
carries on with fresh work. Nothing is lost and nothing wedges -- the
work was redone by somebody else, which is what you asked for -- but
the client's log will show a handful of "File already exists" errors.
This is not new: the same thing happens whenever
resubmit_timed_out_wus() reassigns the workunit of a client that is
merely slow. Reclaiming simply makes it reachable on demand.

Machines leave in groups -- a rack, a reservation, a whole cluster --
so there is a group form of the same action:

    POST /api/v1/clients/reclaim
    {"group_by": "cluster", "group": "alpha", "states": ["stale", "gone"]}

`clients` (a list of ids, matched exactly), `group` with `group_by`,
and `states` may be combined; whatever is given is intersected. A body
that names nothing in particular is refused with 400 rather than taken
to mean every client there is. At most 2000 workunits are marked per
call, and the answer says so when it stopped there.

The marking is done with one `UPDATE ... WHERE status = ASSIGNED AND
wuid IN (...)` per chunk rather than one statement per workunit. Every
`WuAccess.set_status()` takes an EXCLUSIVE lock, and a mass action
taking a thousand of those in a row would stall the thing it is meant
to help. The status is still part of the WHERE, so a workunit that
came back in the meantime is left alone.

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

## Looking at a computation that is not running

    cado-nfs.py --ui-only /tmp/c120/c120.parameters_snapshot.0

serves the monitoring api and the dashboard for an existing working
directory and does nothing else. The read endpoints behave as usual.
Workunit requests are answered 410 and uploads 409, so a client that
happens to still be pointed at that address is told to stop rather than
handed stale work; and the actions are refused with 409, because
setting `NEED_RESUBMIT` with no task running would leave a trap for
whichever task started next.

## Large pools: what a client is, and rolling them up

A pool of any size is unreadable one row at a time, and the useful
unit on a cluster is rarely the individual process: several clients
share a machine, and many machines share a cluster.

    GET /api/v1/clients?group_by=host|cluster|domain
    GET /api/v1/clients?group=<key>&group_by_key=cluster
    GET /api/v1/clients?state=stale,gone

The roll-up is bounded by the number of groups rather than of clients,
which is what makes it usable when the flat list is not. `state`
narrows the list to the clients in the given liveness states; `counts`
and `pool_total` keep describing the whole pool, so "42 of 1400" needs
one request rather than two.

Every answer also carries a `groupings` census:

    "groupings": {"host": 40, "cluster": 3, "domain": 1}

three integers saying how many distinct machines, clusters and domains
the pool covers. The dashboard uses it to pick a roll-up before asking
for one -- the outermost kind that both discriminates and compresses,
since grouping by domain says nothing when there is one domain, and
grouping by machine says nothing when each machine runs one client.

Grouping works with no client cooperation at all, because cado-nfs
names clients predictably: several on one host get `+1`, `+2`
suffixes, and one given no `--clientid` takes `<hostname>.<random>`.
The cluster is then the first dotted component with a trailing
`-<digits>` removed — the same rule `local.sh` uses to name build
trees. `alpha-15+101` is thus machine `alpha-15` in cluster `alpha`.

That is only a guess, though, and it cannot see a domain at all. So a
client introduces itself when it starts:

    POST /clientinfo   {"clientid": ..., "fqdn": ..., "host": ...,
                        "cores": ..., "platform": ...,
                        "overrides": {...}}

which gives the real fully qualified name, the core count, and the
`--override` settings in force — the last of these being how two
clients on one machine that are deliberately configured differently
can be told apart. `self_reported` says which clients did this;
`cado-nfs-client.py` calls it once at startup, best-effort, and an
older server that answers 404 changes nothing.

This endpoint is **unauthenticated**, like the other client endpoints
and for the same reason: a client is given a url and a certificate
fingerprint, never a secret. So what it says is whatever a machine on
`server.whitelist` chose to say. It is used to label and group the
pool for display and for nothing else; fields are coerced, truncated
and whitelisted on the way in, and unrecognised ones are dropped.

If a client reports a hostname, the cluster is derived from *that*
rather than from its client id — saying a machine is `compute-7` but
in cluster `alpha`, because its name happened to start that way,
would be incoherent.

## Drilling into one client or one workunit

    GET /api/v1/clients/<clientid>
    GET /api/v1/workunits/<wuid>

These are the detailed views, and they can afford to be detailed
because they are bounded by construction — one client, one workunit —
unlike the lists, which have a whole pool to worry about.

The client view adds, to what the list already reports, the workunits
it is holding right now, a page of what it has recently handed back
with how long each took, and the raw turnarounds behind its pace
estimate, which the dashboard draws as a sparkline. Whether a machine
is steady, drifting or suddenly much slower says more than any single
number does.

The workunit view adds every attempt at the same piece of work, so a
range that keeps failing shows its whole history and which client had
it each time, and the tail of whatever the client's commands printed —
normally the first thing anybody wants when one has failed. Pass
`output=0` to leave that out.

Attempts are named `<base>` and `<base>__R2`, `__R3` and so on. A
prefix match alone would reach too far, since `c60_sieving_100-200` is
a prefix of `c60_sieving_100-2000`, so the prefix only narrows the
query and the remainder is then checked properly. The captured output
is read off disk, so the path from the database is checked to be
inside the working directory before anything is opened.

`/api/v1/workunits` takes `assigned_to` and `result_from`, and both
match on **substring**. A client id carries a port or an `--override`
suffix, so exact matching would mean knowing the id before being able
to ask about the machine. `assignedclient` is cleared only when a
workunit returns to AVAILABLE, so `assigned_to` also answers "what has
this machine touched". The wildcards `%` and `_` are not escaped:
there is no portable way to say ESCAPE through this database layer,
and a search box where they work is a search box that behaves as one
expects. Reclaiming, by contrast, matches a client id exactly --
reclaiming is not searching.

## Changing a parameter while the computation runs

Almost nothing may be changed. A cado-nfs run is meant to be
reproducible from its `<name>.parameters_snapshot.<N>`, and a knob
turned at three in the morning that appears in no snapshot would break
that quietly, which is the worst way for it to break.

Exactly two parameters are tunable, and they were chosen because
neither can affect what the run computes:

| Parameter | Counter | What it does |
|---|---|---|
| `maxtimedout` | `wu_timedout` | aborts the computation when this many workunits have timed out |
| `maxfailed` | `wu_failed` | aborts the computation when this many have failed |

Both do *nothing at all* except decide when to give up. They change no
workunit, no range and no relation, so raising one cannot alter the
result — only whether the run survives long enough to produce it. That
is also why they are worth exposing: today, when a long run trips one
of these, it dies and you restart it.

    GET  /api/v1/parameters              what they are, and how close
    POST /api/v1/parameters/maxtimedout  {"value": 500}

Everything else — `rels_wanted`, `qrange`, `lpb*`, `admin`/`admax`,
`wutimeout` — is refused with 404, and is meant to be. Note in
particular that `rels_wanted` lives in the task's own state and already
has `request_more_relations()` to mediate it.

Two things keep the record straight:

* **The change is written to the next snapshot in the sequence.** The
  api reads the highest-numbered `parameters_snapshot`, applies the
  override and writes the next one, so "the highest-numbered snapshot
  describes the parameters in force" stays true and resuming from it
  reproduces the configuration rather than the one the run started
  with. If no snapshot can be written the api says so, loudly, in its
  answer and in the log.
* **The change is logged**, by the api when it is made and by each task
  the first time it notices.

A value at or below where the counter already stands is refused with
409: it would abort the computation at the next occurrence, which is
not what anybody means by raising a ceiling.

Tasks read these through `ClientServerTask.tunable()`, which goes to
the database every time rather than through `make_db_dict()` — the
cached flavour would keep serving the value read at startup, which is
precisely the case this exists for.

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
