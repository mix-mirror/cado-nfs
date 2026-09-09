"""
Read-only views of the cado-nfs database, for the monitoring api.

The api server sees the database and nothing else -- not the task
objects, not the scheduler. Everything the monitoring endpoints report
is therefore assembled here out of two sources: the workunits table, and
the DB-backed state dictionaries that the tasks keep (see
cadotask.HasState, and the api_progress table that
cadotask.CompleteFactorization publishes).

Three properties of the database layer shape this module.

 - Attaching to a state dictionary with DictDbDirectAccess takes an
   EXCLUSIVE lock, because its constructor creates the table if it is
   missing. A dashboard polling every couple of seconds must not do
   that, so reads here go through plain SELECTs under a READONLY
   transaction, and a missing table simply reads as empty.

 - The MySQL cursor rewrites some words on their way to the server:
   \\bASC\\b becomes AUTO_INCREMENT, and \\bpurge\\b becomes purgetable.
   Raw SQL below therefore contains neither. Ordering is spelled DESC or
   left implicit.

 - The parameter placeholder is "?" on sqlite and "%s" on MySQL, so it
   is always taken from cursor.parameter_auto_increment rather than
   written out.
"""

import json
import logging
import threading
import time
from datetime import datetime, timezone

from cadofactor import wudb
from cadofactor.database import DictDbDirectAccess
from cadofactor.database import READONLY
from cadofactor.workunit import STATUS_NAMES, WuStatus

logger = logging.getLogger("API server")

# Table in which the api server records whether it is still handing
# workunits out. This lives in the database rather than in the flask
# application object because it is written by cado-nfs.py's main thread
# and read by request handlers, which need not share an address space.
SERVER_STATE_TABLE = "api_server_state"

# Table that cadotask.CompleteFactorization publishes progress into.
PROGRESS_TABLE = "api_progress"

# Statuses that mean "this client returned something".
DONE_STATUSES = (WuStatus.RECEIVED_OK,
                 WuStatus.VERIFIED_OK)
FAILED_STATUSES = (WuStatus.RECEIVED_ERROR,
                   WuStatus.VERIFIED_ERROR)

# How long a cheap aggregate stays good enough to serve again.
AGGREGATE_TTL = 5.0

# Fallback when the computation has not told us its tasks.wutimeout.
DEFAULT_WUTIMEOUT = 10800.0


def status_name(status):
    """
    Name of a numeric workunit status.

    >>> status_name(0)
    'AVAILABLE'
    >>> status_name(5)
    'VERIFIED_OK'
    >>> status_name(99)
    'UNKNOWN(99)'
    """
    try:
        return STATUS_NAMES[int(status)]
    except (IndexError, ValueError, TypeError):
        return "UNKNOWN(%s)" % status


def parse_dbtime(text):
    """
    Parse a timestamp as the workunits table stores it, into epoch
    seconds. The stored form is str(datetime.utcnow()), which drops the
    fractional part when it happens to be zero.

    >>> parse_dbtime('2026-09-09 12:00:00') > 0
    True
    >>> a = parse_dbtime('2026-09-09 12:00:00')
    >>> b = parse_dbtime('2026-09-09 12:00:01.500000')
    >>> round(b - a, 1)
    1.5
    >>> parse_dbtime(None) is None
    True
    >>> parse_dbtime('not a date') is None
    True
    """
    if not text:
        return None
    for fmt in ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S"):
        try:
            naive = datetime.strptime(str(text), fmt)
        except ValueError:
            continue
        # The tasks write utcnow(), so these are UTC despite being naive.
        return naive.replace(tzinfo=timezone.utc).timestamp()
    return None


def liveness(age, in_flight, wutimeout):
    """
    Classify a client from how long ago we last heard from it.

    A client is working if it holds a workunit and we heard from it
    recently; idle if we heard from it recently but it holds nothing;
    stale once it has been silent for longer than the workunit timeout,
    which is the point at which its work starts being reassigned; and
    gone well past that.

    >>> liveness(10, 1, 3600)
    'working'
    >>> liveness(10, 0, 3600)
    'idle'
    >>> liveness(4000, 1, 3600)
    'stale'
    >>> liveness(20000, 1, 3600)
    'gone'
    >>> liveness(None, 0, 3600)
    'unknown'
    """
    if age is None:
        return "unknown"
    if age > 3 * wutimeout:
        return "gone"
    if age > wutimeout:
        return "stale"
    return "working" if in_flight else "idle"


def split_wuid(wuid, name, task_names):
    """
    Split a workunit id into (task, identifier, attempt).

    Workunit ids are built by cadotask.Task.make_wuname as
    <name>_<task>_<identifier>[__R<attempt>]. Both the computation name
    and the task name may themselves contain underscores, so the task is
    recovered by matching against the known task names rather than by
    counting separators.

    >>> tasks = ['polyselect', 'sieve', 'reconstructlog']
    >>> split_wuid('c60_sieve_100-200', 'c60', tasks)
    ('sieve', '100-200', 1)
    >>> split_wuid('c60_sieve_100-200__R3', 'c60', tasks)
    ('sieve', '100-200', 3)
    >>> split_wuid('my_c60_polyselect_0-10', 'my_c60', tasks)
    ('polyselect', '0-10', 1)
    >>> split_wuid('c60_unknowntask_1-2', 'c60', tasks)
    (None, None, 1)
    """
    attempt = 1
    rest = wuid
    marker = rest.rfind("__R")
    if marker >= 0:
        tail = rest[marker + 3:]
        if tail.isdigit():
            attempt = int(tail)
            rest = rest[:marker]
    if name and rest.startswith(name + "_"):
        rest = rest[len(name) + 1:]
    for task in sorted(task_names, key=len, reverse=True):
        if rest == task:
            return (task, "", attempt)
        if rest.startswith(task + "_"):
            return (task, rest[len(task) + 1:], attempt)
    return (None, None, attempt)


class _PerThread(object):
    """
    Keeps one object per thread.

    sqlite3 connections may only be used from the thread that created
    them, which is why ApiServer keeps a per-thread connection pool.
    Anything built on top of such a connection inherits the constraint.
    """

    def __init__(self, factory):
        self._factory = factory
        self._objects = {}
        self._lock = threading.Lock()

    def get(self):
        tid = threading.current_thread().ident
        with self._lock:
            obj = self._objects.get(tid)
            if obj is None:
                obj = self._factory()
                self._objects[tid] = obj
            return obj


class ServingState(object):
    """
    The "are we still handing workunits out" flag, kept in the database.

    cado-nfs.py's main thread clears it when the distributed phase is
    over; request handlers read it when a client asks for work. The two
    are in different threads today and would be in different processes
    under gunicorn, so an attribute on the flask application will not do.

    Reads are cached briefly, because they sit on the hot path of every
    single workunit request.
    """

    KEY = "serving_wus"
    TTL = 1.0

    def __init__(self, connection_factory):
        self._dicts = _PerThread(
            lambda: DictDbDirectAccess(connection_factory(),
                                       SERVER_STATE_TABLE))
        self._cached = None
        self._cached_at = 0.0

    def set(self, serving):
        self._dicts.get()[self.KEY] = bool(serving)
        self._cached = bool(serving)
        self._cached_at = time.time()

    def get(self):
        now = time.time()
        if self._cached is not None and now - self._cached_at < self.TTL:
            return self._cached
        try:
            value = self._dicts.get().get(self.KEY, True)
        except Exception as e:
            logger.warning("Could not read serving state (%s),"
                           " assuming we still serve workunits", e)
            value = True
        self._cached = bool(value)
        self._cached_at = now
        return self._cached


class DbViews(object):
    """
    The read-only queries that the monitoring endpoints are built from.

    Every method here takes a connection from the factory it was given,
    so that it is usable from any request thread.
    """

    def __init__(self, connection_factory):
        self._connection_factory = connection_factory
        self._cache = {}
        self._lock = threading.Lock()

    # ---------------- plumbing ----------------

    def _read(self, function, *args, **kwargs):
        """
        Run function(cursor, ...) under a read-only transaction.
        """
        connection = self._connection_factory()
        return connection.harness_transaction(READONLY, function,
                                              *args, **kwargs)

    def _cached(self, key, producer, ttl=AGGREGATE_TTL):
        """
        Memoize an aggregate for a short while.

        The aggregates below are a group-by over the whole workunits
        table. They are cheap enough with the indices that wudb.WuTable
        declares, but not so cheap that a dashboard should recompute
        them for every panel it draws.
        """
        now = time.time()
        with self._lock:
            entry = self._cache.get(key)
            if entry is not None and now - entry[0] < ttl:
                return entry[1]
        value = producer()
        with self._lock:
            self._cache[key] = (now, value)
        return value

    def invalidate(self):
        """
        Drop memoized aggregates, after an action that changed them.
        """
        with self._lock:
            self._cache.clear()

    # ---------------- state dictionaries ----------------

    def read_state_table(self, name):
        """
        Read a DB-backed state dictionary by plain SELECT.

        Unlike DictDbDirectAccess this neither creates the table nor
        takes an exclusive lock; a table that does not exist reads as an
        empty dictionary, which is the right answer for a task that has
        not run yet.
        """
        try:
            wudb.check_tablename(name)
        except Exception:
            logger.warning("refusing to read bogus table name %r", name)
            return {}

        def query(cursor):
            cursor.execute("SELECT kkey, type, value FROM %s;" % name)
            return cursor.cursor.fetchall()

        try:
            rows = self._read(query)
        except Exception as e:
            # Most often: the task never ran, so its table is absent.
            logger.debug("could not read state table %s (%s)", name, e)
            return {}

        types = DictDbDirectAccess.types
        out = {}
        for key, typeindex, value in rows:
            try:
                constructor = types[int(typeindex)]
                if constructor is bool:
                    out[key] = (value == "True")
                else:
                    out[key] = constructor(value)
            except (IndexError, ValueError, TypeError):
                out[key] = value
        return out

    def progress_state(self):
        """
        The api_progress table, with its JSON fields decoded.
        """
        raw = self.read_state_table(PROGRESS_TABLE)
        out = dict(raw)
        for key in ("pipeline", "done"):
            try:
                out[key] = json.loads(raw.get(key, "[]"))
            except ValueError:
                out[key] = []
        return out

    def task_names(self):
        return [t.get("name") for t in self.progress_state().get(
            "pipeline", []) if t.get("name")]

    def wutimeout(self):
        try:
            return float(self.progress_state().get("wutimeout")
                         or DEFAULT_WUTIMEOUT)
        except (TypeError, ValueError):
            return DEFAULT_WUTIMEOUT

    # ---------------- workunits ----------------

    def status_counts(self):
        """
        Number of workunits in each status, keyed by status name.
        """
        def produce():
            def query(cursor):
                cursor.execute("SELECT status, COUNT(*) FROM workunits"
                               " GROUP BY status;")
                return cursor.cursor.fetchall()

            counts = {name: 0 for name in STATUS_NAMES}
            for status, count in self._read(query):
                counts[status_name(status)] = int(count)
            return counts

        return self._cached("status_counts", produce)

    def client_aggregates(self):
        """
        Per-client tallies, as two group-by queries.

        The first covers what is in flight right now, and seeks on
        status; the second covers what each client has returned, and is
        answered from the (resultclient, status) index. Neither reads
        the workunit bodies.
        """
        def produce():
            def in_flight(cursor):
                qm = cursor.parameter_auto_increment
                cursor.execute(
                    "SELECT assignedclient, COUNT(*),"
                    " MIN(timeassigned), MAX(timeassigned)"
                    " FROM workunits WHERE status = " + qm +
                    " GROUP BY assignedclient;",
                    [int(WuStatus.ASSIGNED)])
                return cursor.cursor.fetchall()

            def returned(cursor):
                cursor.execute(
                    "SELECT resultclient, status, COUNT(*),"
                    " MAX(timeresult)"
                    " FROM workunits WHERE resultclient IS NOT NULL"
                    " GROUP BY resultclient, status;")
                return cursor.cursor.fetchall()

            clients = {}

            def entry(name):
                return clients.setdefault(name, {
                    "clientid": name,
                    "in_flight": 0,
                    "completed": 0,
                    "failed": 0,
                    "oldest_assignment": None,
                    "last_assignment": None,
                    "last_result": None,
                })

            for name, count, oldest, newest in self._read(in_flight):
                if name is None:
                    continue
                e = entry(name)
                e["in_flight"] = int(count)
                e["oldest_assignment"] = parse_dbtime(oldest)
                e["last_assignment"] = parse_dbtime(newest)

            for name, status, count, newest in self._read(returned):
                if name is None:
                    continue
                e = entry(name)
                if int(status) in DONE_STATUSES:
                    e["completed"] += int(count)
                elif int(status) in FAILED_STATUSES:
                    e["failed"] += int(count)
                stamp = parse_dbtime(newest)
                if stamp is not None and (e["last_result"] is None
                                          or stamp > e["last_result"]):
                    e["last_result"] = stamp
            return clients

        return self._cached("clients", produce)

    def clients(self):
        """
        Per-client view, with liveness and contribution share.
        """
        aggregates = self.client_aggregates()
        timeout = self.wutimeout()
        now = time.time()
        total_completed = sum(c["completed"] for c in aggregates.values())

        out = []
        for entry in aggregates.values():
            client = dict(entry)
            stamps = [s for s in (entry["last_assignment"],
                                  entry["last_result"]) if s is not None]
            last_seen = max(stamps) if stamps else None
            client["last_seen"] = last_seen
            # Only the liveness verdict is computed here, because it
            # needs tasks.wutimeout, which the client does not have. The
            # ages themselves stay out of the body so that an unchanged
            # answer revalidates as a 304; consumers subtract last_seen
            # from the X-Cado-Server-Time header instead.
            age = None if last_seen is None else max(0.0, now - last_seen)
            client["state"] = liveness(age, entry["in_flight"], timeout)
            client["share"] = (entry["completed"] / total_completed
                               if total_completed else 0.0)
            out.append(client)
        out.sort(key=lambda c: (-c["completed"], c["clientid"]))
        return out

    def list_workunits(self, status=None, assigned_to=None,
                       result_from=None, task=None,
                       assigned_older_than=None, limit=50, offset=0):
        """
        A page of the workunits table, most recent first.

        The two client filters are deliberately separate rather than one
        "client" filter: a workunit has both an assignedclient and a
        resultclient, and which one you mean is the difference between
        "what is this machine chewing on" and "what did it hand back".

        Filtering by task uses a LIKE on the wuid prefix, since the task
        name is encoded there and there is no column for it.
        """
        conditions = {}
        equalities = {}
        if status is not None:
            equalities["status"] = int(status)
        if assigned_to is not None:
            equalities["assignedclient"] = assigned_to
        if result_from is not None:
            equalities["resultclient"] = result_from
        if equalities:
            conditions["eq"] = equalities
        if task is not None:
            name = self.progress_state().get("name", "")
            prefix = ("%s_%s_" % (name, task)) if name else ("%s_" % task)
            conditions["like"] = {"wuid": prefix + "%"}
        if assigned_older_than is not None:
            stamp = datetime.fromtimestamp(
                time.time() - float(assigned_older_than),
                timezone.utc).replace(tzinfo=None)
            conditions["lt"] = {"timeassigned": str(stamp)}

        def query(cursor):
            return cursor.where_as_dict("workunits",
                                        limit=int(limit),
                                        offset=int(offset),
                                        order=("wurowid", "DESC"),
                                        **conditions)

        rows = self._read(query)
        return [self.describe_workunit(row) for row in rows]

    def describe_workunit(self, row, with_body=False):
        """
        Turn a workunits row into what the api hands out.
        """
        name = self.progress_state().get("name", "")
        task, identifier, attempt = split_wuid(row["wuid"], name,
                                               self.task_names())
        out = {
            "wuid": row["wuid"],
            "status": int(row["status"]),
            "status_name": status_name(row["status"]),
            "task": task,
            "identifier": identifier,
            "attempt": attempt,
            "assignedclient": row.get("assignedclient"),
            "resultclient": row.get("resultclient"),
            "timecreated": parse_dbtime(row.get("timecreated")),
            "timeassigned": parse_dbtime(row.get("timeassigned")),
            "timeresult": parse_dbtime(row.get("timeresult")),
            "timeverified": parse_dbtime(row.get("timeverified")),
            "errorcode": row.get("errorcode"),
            "failedcommand": row.get("failedcommand"),
        }
        if with_body:
            try:
                out["workunit"] = json.loads(row["wu"])
            except (ValueError, KeyError, TypeError):
                out["workunit"] = None
            files = row.get("files") or []
            out["files"] = [{"filename": f.get("filename"),
                             "path": f.get("path"),
                             "type": f.get("type")} for f in files]
        return out

    # ---------------- assembled views ----------------

    def task_view(self, entry, current, done_names):
        """
        Describe one task of the pipeline from its state dictionary.
        """
        name = entry.get("name")
        state = self.read_state_table(name) if name else {}
        if name == current:
            phase = "running"
        elif name in done_names:
            phase = "done"
        else:
            phase = "pending"

        view = {"name": name,
                "title": entry.get("title"),
                "phase": phase}

        for key in ("achievement", "eta", "progress_time",
                    "wu_submitted", "wu_received", "wu_timedout",
                    "wu_failed", "wu_range_received"):
            if key in state:
                view[key] = state[key]

        # Whatever the task happens to publish about its own subject
        # matter. Absent keys simply mean "not this kind of task".
        highlights = {}
        for key in ("rels_found", "rels_wanted", "qnext", "adnext",
                    "noutrels", "nr_poly_submitted"):
            if key in state:
                highlights[key] = state[key]
        if highlights:
            view["highlights"] = highlights

        times = {k: v for k, v in state.items()
                 if k.startswith("cputime_") or k.startswith("realtime_")}
        if times:
            view["times"] = times
        return view

    def progress(self):
        """
        The whole "where do we stand" answer.
        """
        published = self.progress_state()
        current = published.get("current") or None
        done = published.get("done", [])
        done_names = {d.get("name") for d in done}
        stats = {d.get("name"): d.get("stats", []) for d in done}

        tasks = []
        for entry in published.get("pipeline", []):
            view = self.task_view(entry, current, done_names)
            if view["name"] in stats:
                view["stats"] = stats[view["name"]]
            tasks.append(view)

        started = published.get("current_started") or None
        return {
            "computation": published.get("computation"),
            "algorithm": published.get("algo"),
            "current": current,
            "current_started": started,
            "finished": bool(published.get("finished", False)),
            "elapsed": published.get("elapsed"),
            "cputotal": published.get("cputotal"),
            "tasks": tasks,
        }


if __name__ == "__main__":
    import doctest
    doctest.testmod()
