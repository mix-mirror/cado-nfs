"""
A synthetic cado-nfs database that looks like a computation in
mid-sieving, and an ApiServer bound to it.

Standing up a real factorization just to check that the monitoring api
reports the right numbers would be slow and flaky, and would need the
compiled binaries. Everything the api reads is in the database, so we
write the database directly and assert against that.

The population below is chosen to exercise the awkward cases:

 - three clients in three different states: one busy, one busy but
   quiet for a while, one that has been silent long enough to count as
   gone;
 - a workunit still assigned from a task that has already finished,
   which the reclaim actions must refuse to touch, because
   ClientServerTask.resubmit_timed_out_wus() would charge it to the
   task that is running now;
 - both error statuses, so that the per-client failure tallies have
   something to count.
"""

import json
import os
import time
from datetime import datetime, timedelta

from cadofactor import wudb
from cadofactor.database import EXCLUSIVE, DictDbDirectAccess
from cadofactor.workunit import Workunit, WuStatus
from cadofactor.api_server import ApiServer

NAME = "c60"

# tasks.wutimeout as the fixture publishes it. Client liveness is
# expressed in multiples of this.
WUTIMEOUT = 3600
# Deliberately large. This sets the small-sample staleness floor, and
# therefore how big BRISK's threshold may be while still being *below*
# it -- which is what keeps every client's state far away from a
# boundary. See MIN_STATE_MARGIN.
WUTIMEOUTCHECK = 600

PIPELINE = [
    {"name": "polyselect",
     "title": "Polynomial Selection (size optimized)"},
    {"name": "polyselect2",
     "title": "Polynomial Selection (root optimized)"},
    {"name": "factorbase", "title": "Generate Factor Base"},
    {"name": "freerel", "title": "Generate Free Relations"},
    {"name": "sieving", "title": "Lattice Sieving"},
    {"name": "duplicates1", "title": "Filtering - Duplicate Removal"},
    {"name": "purge", "title": "Filtering - Singleton removal"},
    {"name": "merge", "title": "Filtering - Merging"},
    # tasks.linalg.run=false. A disabled task does not get stepped
    # over: Task.run() raises, so the run stops there and sqrt never
    # happens either. The api is expected to say both things.
    {"name": "linalg", "title": "Linear Algebra", "run": False},
    {"name": "sqrt", "title": "Square Root"},
]

# How the workunits table is filled. Ages are in seconds before now.
#   (wuid suffix, status, assigned_age, client, result_age, resultclient)
BUSY = "alpha-01"
SLOW = "alpha-02"
GONE = "alpha-03"
# A brisk client that has just gone quiet. It is the case that
# distinguishes the per-client staleness threshold from the global one:
# 20 minutes of silence is nothing next to a one-hour wutimeout, but it
# is a long time for a machine that returns a workunit every two
# minutes.
FAST = "alpha-04"
FAST_TURNAROUND = 120
FAST_SILENT_FOR = 1200
# A client so brisk that six times its pace is below the small-sample
# floor. With enough samples agreeing, the floor is dropped and the
# threshold is purely what this client has shown us -- which is the
# case that tells us the sample count is really reaching stale_after.
# The one workunit FAST is holding, which is what a reclaim of the
# stale members of its cluster has to take back.
FAST_WORKUNIT = NAME + "_sieving_650000-651000"
BRISK = "alpha-05"
BRISK_TURNAROUND = 150
BRISK_SAMPLES = 12

# A piece of work that failed once and was retried. RETRY_DECOY exists
# to catch a prefix match that reaches too far: its id begins with
# RETRY_BASE, but it is a different workunit, not another attempt.
RETRY_BASE = NAME + "_sieving_990000-991000"
RETRY_SECOND = RETRY_BASE + "__R2"
RETRY_DECOY = NAME + "_sieving_990000-9910000"

# What the tests expect to find, so that a change to the population
# above shows up as a failure here rather than as a silent drift.
# No client's state may be closer than this to flipping. The fixture
# describes ages relative to "now", so a slow or loaded machine will
# read it later than it was written; a test that only holds for a few
# seconds is a test that passes on the author's laptop and fails on
# somebody else's runner. This must comfortably exceed the ctest
# timeout, since a run that takes longer than that is killed anyway.
MIN_STATE_MARGIN = 300

EXPECT = {
    "current_task": "sieving",
    "clients": {
        BUSY: {"state": "working", "in_flight": 4, "completed": 92,
               "failed": 0},
        SLOW: {"state": "working", "in_flight": 1, "completed": 60,
               "failed": 3},
        GONE: {"state": "gone", "in_flight": 4, "completed": 28,
               "failed": 0},
        FAST: {"state": "stale", "in_flight": 1, "completed": 12,
               "failed": 0, "liveness_basis": "turnaround"},
        BRISK: {"state": "idle", "in_flight": 0,
                "completed": BRISK_SAMPLES, "failed": 0,
                "liveness_basis": "turnaround"},
    },
    # What the old, global-only rule would have said about FAST. The
    # point of the per-client threshold is that these differ.
    "fast_state_under_wutimeout_only": "working",
    # Reclaiming GONE must take its three sieving workunits and leave
    # the polyselect leftover alone.
    "reclaimable_from_gone": 3,
    "stale_from_gone": NAME + "_polyselect_5000-5100",
}


def state_margin(client, now):
    """
    How long a client may go unread before its reported state changes.

    liveness() switches at `stale_after` and at three times it, so the
    margin is the distance from the client's current age to whichever
    of those boundaries is next -- or, for a client already past the
    last one, infinity. Ages are not in the api's answers by design, so
    they are recomputed here from last_seen.
    """
    import math
    if client.get("last_seen") is None:
        return math.inf
    age = now - client["last_seen"]
    threshold = client["stale_after"]
    for boundary in (threshold, 3 * threshold):
        if age < boundary:
            return boundary - age
    return math.inf


def utc(age):
    return str(datetime.utcnow() - timedelta(seconds=age))


def populate(db, workdir):
    wuaccess = wudb.WuAccess(db)
    wuaccess.create_tables()
    conn = db.connect()

    tasks = DictDbDirectAccess(conn, "tasks")
    tasks.update({"workdir": workdir,
                  "N": 90377629292003121684002147101760858109247336549,
                  "starttime": time.time() - 7200.0,
                  "elapsed": 0.0})

    progress = DictDbDirectAccess(conn, "api_progress")
    progress.update({
        "pipeline": json.dumps(PIPELINE),
        "computation": "factorization",
        "algo": "nfs",
        "name": NAME,
        "wutimeout": WUTIMEOUT,
        "wutimeoutcheck": WUTIMEOUTCHECK,
        "tunable_defaults": json.dumps({"maxtimedout": 100,
                                        "maxfailed": 100}),
        "current": "sieving",
        "current_started": time.time() - 5400.0,
        "done": json.dumps([
            {"name": "polyselect", "time": time.time() - 6000,
             "stats": ["Total cpu/real time for polyselect: 812.3/205.1"]},
            {"name": "polyselect2", "time": time.time() - 5800},
            {"name": "factorbase", "time": time.time() - 5600},
            {"name": "freerel", "time": time.time() - 5500},
        ]),
        "finished": False})

    DictDbDirectAccess(conn, "polyselect").update({
        "wu_submitted": 40, "wu_received": 40, "wu_timedout": 0,
        "wu_failed": 0, "wu_range_received": 4000,
        "achievement": 1.0, "cputime_polyselect": 812.3})

    DictDbDirectAccess(conn, "sieving").update({
        "wu_submitted": 220, "wu_received": 181, "wu_timedout": 3,
        "wu_failed": 2, "wu_range_received": 18100,
        "rels_found": 1943221, "rels_wanted": 3110000,
        "qnext": 1810000,
        "achievement": 0.6248,
        "eta": "Wed Sep  9 18:42:11 2026",
        "progress_time": time.time() - 45,
        "cputime_las": 194322.5, "realtime_las": 24122.0})

    rows = []

    def add(wuid, status, assigned=None, client=None,
            result_age=None, resultclient=None, errorcode=None):
        rows.append((wuid, status, assigned, client, result_age,
                     resultclient, errorcode))

    for i in range(12):
        add("%s_sieving_%d-%d" % (NAME, 900000 + i * 1000,
                                  901000 + i * 1000), WuStatus.AVAILABLE)

    for i in range(4):
        add("%s_sieving_%d-%d" % (NAME, 950000 + i * 1000,
                                  951000 + i * 1000),
            WuStatus.ASSIGNED, assigned=60 + i * 30, client=BUSY)

    add("%s_sieving_970000-971000" % NAME, WuStatus.ASSIGNED,
        assigned=300, client=SLOW)

    # Silent for longer than 3 * wutimeout, so: gone.
    for i in range(3):
        add("%s_sieving_%d-%d" % (NAME, 980000 + i * 1000,
                                  981000 + i * 1000),
            WuStatus.ASSIGNED, assigned=4 * WUTIMEOUT + i, client=GONE)

    # The leftover from a finished task. Reclaiming must skip it.
    add(EXPECT["stale_from_gone"], WuStatus.ASSIGNED,
        assigned=5 * WUTIMEOUT, client=GONE)

    # The age of the *most recent* result matters: liveness is decided
    # on the last time we heard from a client at all, so the one that
    # is supposed to look gone must not have handed anything in
    # recently either.
    for count, client, base, newest in ((90, BUSY, 100000, 3000),
                                        (60, SLOW, 400000, 3000),
                                        (28, GONE, 700000,
                                         5 * WUTIMEOUT)):
        for i in range(count):
            add("%s_sieving_%d-%d" % (NAME, base + i * 1000,
                                      base + 1000 + i * 1000),
                WuStatus.VERIFIED_OK,
                assigned=newest + 1000 + i * 20, client=client,
                result_age=newest + i * 20, resultclient=client)

    # FAST: a short, very regular turnaround, and then silence.
    for i in range(12):
        add("%s_sieving_%d-%d" % (NAME, 600000 + i * 1000,
                                  601000 + i * 1000),
            WuStatus.VERIFIED_OK,
            assigned=FAST_SILENT_FOR + FAST_TURNAROUND + i * 300,
            client=FAST,
            result_age=FAST_SILENT_FOR + i * 300, resultclient=FAST)
    add(FAST_WORKUNIT, WuStatus.ASSIGNED,
        assigned=FAST_SILENT_FOR, client=FAST)

    for i in range(BRISK_SAMPLES):
        add("%s_sieving_%d-%d" % (NAME, 300000 + i * 1000,
                                  301000 + i * 1000),
            WuStatus.VERIFIED_OK,
            assigned=4 + BRISK_TURNAROUND + i * 10, client=BRISK,
            result_age=4 + i * 10, resultclient=BRISK)

    # One workunit that failed and was retried, plus the decoy whose
    # id merely starts with the same text.
    add(RETRY_BASE, WuStatus.VERIFIED_ERROR, assigned=900, client=SLOW,
        result_age=850, resultclient=SLOW, errorcode=1)
    add(RETRY_SECOND, WuStatus.VERIFIED_OK, assigned=800, client=BUSY,
        result_age=700, resultclient=BUSY)
    add(RETRY_DECOY, WuStatus.VERIFIED_OK, assigned=600, client=BUSY,
        result_age=560, resultclient=BUSY)

    for i in range(2):
        add("%s_sieving_%d-%d" % (NAME, 800000 + i * 1000,
                                  801000 + i * 1000),
            WuStatus.VERIFIED_ERROR, assigned=5000, client=SLOW,
            result_age=4900, resultclient=SLOW, errorcode=1)

    table = wuaccess.mapper.table

    def insert(cursor):
        for (wuid, status, assigned, client, result_age, resultclient,
             errorcode) in rows:
            d = {"wuid": wuid,
                 "wu": str(Workunit(id=wuid, commands=["true"],
                                    timeout=3600)),
                 "status": status,
                 "timecreated": utc(20000)}
            if assigned is not None:
                d["timeassigned"] = utc(assigned)
                d["assignedclient"] = client
            if result_age is not None:
                d["timeresult"] = utc(result_age)
                d["resultclient"] = resultclient
                d["timeverified"] = utc(result_age - 1)
            if errorcode is not None:
                d["errorcode"] = errorcode
            table.insert(cursor, d)

    conn.harness_transaction(EXCLUSIVE, insert)

    # A captured stderr for the attempt that failed, so that the detail
    # view has something to show. This is what the server itself writes
    # when a client uploads.
    upload = os.path.join(workdir, NAME + ".upload")
    if not os.path.isdir(upload):
        os.makedirs(upload)
    stderr_name = RETRY_BASE + ".stderr0"
    stderr_path = os.path.join(upload, stderr_name)
    with open(stderr_path, "w") as f:
        for i in range(200):
            f.write("las: some chatter on line %d\n" % i)
        f.write("las: Error: sieve region is empty, giving up\n")

    def attach(cursor):
        wuaccess._add_files(cursor,
                            [(stderr_name, stderr_path, "stderr0", 0)],
                            wuid=RETRY_BASE)

    conn.harness_transaction(EXCLUSIVE, attach)
    del tasks, progress


def write_snapshot(workdir):
    """
    A parameters snapshot for the api to build the next one on.
    """
    path = os.path.join(workdir, NAME + ".parameters_snapshot.0")
    with open(path, "w") as f:
        f.write("name = %s\n" % NAME)
        f.write("tasks.workdir = %s\n" % workdir)
        f.write("tasks.wutimeout = %d\n" % WUTIMEOUT)
    return path


def write_log(workdir):
    path = os.path.join(workdir, NAME + ".log")
    with open(path, "w") as f:
        for i in range(300):
            f.write("PID1 2026-09-09 12:00:%02d,000 Info:Lattice Sieving:"
                    " line %d\n" % (i % 60, i))
        f.write("PID1 2026-09-09 12:05:00,000 Error:Lattice Sieving:"
                " something went wrong\n")
    return path


def build(workdir, **kwargs):
    """
    Create the database and return an ApiServer bound to it.
    """
    if not os.path.isdir(workdir):
        os.makedirs(workdir)
    db = wudb.DBFactory("db:sqlite3://%s/%s.db" % (workdir, NAME),
                        create=True)
    populate(db, workdir)
    write_log(workdir)
    write_snapshot(workdir)
    options = dict(threaded=False,
                   uploaddir=os.path.join(workdir, NAME + ".upload"),
                   nrsubdir=0,
                   cafile=None,
                   whitelist=["127.0.0.1/32"],
                   workdir=workdir,
                   name=NAME)
    options.update(kwargs)
    return ApiServer(None, 0, db, **options)
