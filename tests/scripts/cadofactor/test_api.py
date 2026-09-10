#!/usr/bin/env python3

"""
Checks on the cado-nfs monitoring api, against the synthetic database
that api_fixture.py builds.

Run as:  test_api.py <section> <workdir>
where section is one of openapi, auth, views, actions, or serve (which
starts a real http server for test_monitor_cli.sh to talk to).
"""

import os
import stat
import sys
import threading

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import api_fixture                                          # noqa: E402
from api_fixture import (EXPECT, NAME, WUTIMEOUT, FAST,     # noqa: E402
                         FAST_TURNAROUND, FAST_SILENT_FOR)  # noqa: E402


class Failures(object):
    def __init__(self):
        self.count = 0

    def check(self, condition, description, detail=None):
        if condition:
            print("ok       %s" % description)
        else:
            self.count += 1
            print("FAILED   %s" % description)
            if detail is not None:
                print("         %s" % (detail,))

    def equal(self, got, want, description):
        self.check(got == want, description,
                   None if got == want else "got %r, want %r" % (got, want))


def authorized(app):
    return {"Authorization": "Bearer " + app.api_token}


# ------------------------------------------------------------------
# openapi
# ------------------------------------------------------------------


def section_openapi(app, f):
    from cadofactor.api.spec import rule_to_openapi_path

    client = app.test_client()
    response = client.get("/api/v1/openapi.json")
    f.equal(response.status_code, 200,
            "the OpenAPI document is served without a token")
    doc = response.get_json()

    f.equal(doc.get("openapi"), "3.1.0", "it declares OpenAPI 3.1.0")
    f.check("info" in doc and "title" in doc["info"], "it has an info title")
    f.check(doc.get("components", {}).get("securitySchemes", {})
            .get("bearerAuth", {}).get("scheme") == "bearer",
            "it declares the bearer security scheme")

    # Every route that is reachable must be described. This is what
    # keeps the document from drifting: routes and documentation come
    # from the same @api_route declaration, and anything registered by
    # hand behind its back shows up here.
    reachable = set()
    for rule in app.url_map.iter_rules():
        if rule.endpoint == "static":
            continue
        reachable.add(rule_to_openapi_path(str(rule)))
    documented = set(doc["paths"])
    f.equal(sorted(reachable - documented), [],
            "every reachable route appears in the document")
    f.equal(sorted(documented - reachable), [],
            "every documented path is reachable")

    # And every operation must actually say something.
    incomplete = []
    for path, entry in doc["paths"].items():
        for method, op in entry.items():
            if not op.get("summary"):
                incomplete.append("%s %s: no summary" % (method, path))
            if not op.get("operationId"):
                incomplete.append("%s %s: no operationId" % (method, path))
            if not op.get("responses"):
                incomplete.append("%s %s: no responses" % (method, path))
            if not op.get("tags"):
                incomplete.append("%s %s: no tags" % (method, path))
    f.equal(incomplete, [], "every operation is described")

    ids = [op["operationId"] for e in doc["paths"].values()
           for op in e.values()]
    f.equal(len(ids), len(set(ids)), "operationIds are unique")

    # The endpoints that need a token must say so, and the ones clients
    # use must not.
    secured = {(p, m) for p, e in doc["paths"].items()
               for m, op in e.items() if op.get("security")}
    f.check(("/api/v1/clients", "get") in secured,
            "monitoring endpoints are marked as needing a token")
    f.check(("/workunit", "get") not in secured,
            "the workunit endpoint is not marked as needing a token")


# ------------------------------------------------------------------
# auth
# ------------------------------------------------------------------


def section_auth(app, f):
    client = app.test_client()

    f.equal(client.get("/api/v1/info").status_code, 401,
            "no token is refused")
    f.equal(client.get("/api/v1/info",
                       headers={"Authorization": "Bearer wrong"}
                       ).status_code, 401,
            "a wrong token is refused")
    response = client.get("/api/v1/info", headers=authorized(app))
    f.equal(response.status_code, 200, "the right token is accepted")

    f.check("Bearer" in client.get("/api/v1/info").headers.get(
        "WWW-Authenticate", ""),
        "a refusal carries a WWW-Authenticate header")

    # Client endpoints must keep working with no token at all: clients
    # are never given one.
    f.equal(client.get("/").status_code, 200,
            "the liveness probe needs no token")
    f.equal(client.get("/files").status_code, 200,
            "the file list needs no token")
    f.equal(client.get("/workunit", data={"clientid": "t"}).status_code,
            200, "fetching a workunit needs no token")

    if os.name == "posix":
        mode = stat.S_IMODE(os.stat(app.token_file).st_mode)
        f.equal(oct(mode), oct(0o600), "the token file is mode 0600")

        # A token anybody can read is not a token.
        from cadofactor.api import auth
        os.chmod(app.token_file, 0o644)
        try:
            auth.load_or_create_token(app.token_file)
            f.check(False, "a world-readable token file is refused")
        except auth.TokenError:
            f.check(True, "a world-readable token file is refused")
        finally:
            os.chmod(app.token_file, 0o600)

        # And an existing private one is reused, so that a monitor
        # survives a restart of cado-nfs.py.
        f.equal(auth.load_or_create_token(app.token_file), app.api_token,
                "an existing private token file is reused")

    # An error carries the status the exception says, not a blanket 404:
    # cado-nfs-client.py keys its clean shutdown off the 410 below.
    app.serving.set(False)
    f.equal(client.get("/workunit", data={"clientid": "t"}).status_code,
            410, "a client is told 410 once we stop serving workunits")
    app.serving.set(True)
    f.equal(client.get("/file/not-registered").status_code, 404,
            "an unregistered file is 404, not 500")


# ------------------------------------------------------------------
# views
# ------------------------------------------------------------------


def section_views(app, f):
    client = app.test_client()
    headers = authorized(app)

    def get(path):
        response = client.get(path, headers=headers)
        f.equal(response.status_code, 200, "GET %s" % path)
        return response

    info = get("/api/v1/info").get_json()
    f.equal(info["name"], NAME, "info reports the computation name")
    f.equal(info["computation"], "factorization",
            "info reports the kind of computation")
    f.equal(info["wutimeout"], float(WUTIMEOUT),
            "info reports tasks.wutimeout")

    progress = get("/api/v1/progress").get_json()
    f.equal(progress["current"], EXPECT["current_task"],
            "progress names the running task")
    phases = {t["name"]: t["phase"] for t in progress["tasks"]}
    f.equal(phases.get("sieving"), "running", "sieving is running")
    f.equal(phases.get("polyselect"), "done", "polyselect is done")
    f.equal(phases.get("merge"), "pending", "merge is pending")
    f.equal(len(progress["tasks"]), len(api_fixture.PIPELINE),
            "the whole pipeline is reported")

    sieving = [t for t in progress["tasks"] if t["name"] == "sieving"][0]
    f.equal(round(sieving["achievement"], 4), 0.6248,
            "the task's own achievement is passed through")
    f.equal(sieving["highlights"]["rels_found"], 1943221,
            "sieving reports the relations it has found")
    f.equal(sieving["highlights"]["rels_wanted"], 3110000,
            "sieving reports the relations it wants")
    f.check("eta" in sieving, "sieving reports an ETA")

    summary = get("/api/v1/workunits/summary").get_json()
    counts = summary["counts"]
    f.equal(counts["AVAILABLE"], 12, "available workunits are counted")
    f.equal(counts["ASSIGNED"], 10, "assigned workunits are counted")
    f.equal(counts["VERIFIED_OK"], 188, "finished workunits are counted")
    f.equal(counts["VERIFIED_ERROR"], 2, "failed workunits are counted")
    f.equal(summary["total"], sum(counts.values()), "the total adds up")

    clients = {c["clientid"]: c
               for c in get("/api/v1/clients").get_json()["clients"]}
    f.equal(sorted(clients), sorted(EXPECT["clients"]),
            "every client that has been seen is reported")
    for name, want in EXPECT["clients"].items():
        got = clients.get(name, {})
        for key, value in want.items():
            f.equal(got.get(key), value, "%s: %s" % (name, key))
    f.check(abs(sum(c["share"] for c in clients.values()) - 1.0) < 1e-9,
            "the contribution shares add up to one")

    # The staleness threshold is per client, learnt from how long its
    # workunits have recently been taking. FAST is the case that
    # separates that from the old global rule: it has been silent for
    # 20 minutes, which is nothing against a one-hour wutimeout but a
    # long time for a machine that returns a workunit every two.
    from cadofactor.api import views
    fast = clients[FAST]
    f.check(fast["turnaround_samples"] >= views.TURNAROUND_MIN_SAMPLES,
            "%s: enough samples to judge its turnaround" % FAST)
    f.check(abs(fast["typical_turnaround"] - FAST_TURNAROUND) < 1.0,
            "%s: its typical turnaround is measured" % FAST,
            "got %r, want about %d" % (fast["typical_turnaround"],
                                       FAST_TURNAROUND))
    f.equal(fast["liveness_basis"], "turnaround",
            "%s: judged on its own turnaround, not on wutimeout" % FAST)
    f.check(fast["stale_after"] < WUTIMEOUT,
            "%s: is given less rope than wutimeout would" % FAST,
            "stale_after=%r wutimeout=%r" % (fast["stale_after"],
                                             WUTIMEOUT))
    f.equal(views.liveness(FAST_SILENT_FOR, 1, WUTIMEOUT),
            EXPECT["fast_state_under_wutimeout_only"],
            "%s: the global rule alone would have called it %s"
            % (FAST, EXPECT["fast_state_under_wutimeout_only"]))

    # A client we have learnt nothing about falls back to wutimeout,
    # and the threshold never exceeds it.
    for name, entry in clients.items():
        f.check(entry["stale_after"] <= WUTIMEOUT,
                "%s: threshold never exceeds wutimeout" % name)
        if entry["turnaround_samples"] < views.TURNAROUND_MIN_SAMPLES:
            f.equal(entry["liveness_basis"], "wutimeout",
                    "%s: too few samples, so falls back" % name)

    # Ages must not be in the body, or no answer would ever revalidate.
    f.check("idle_seconds" not in clients[list(clients)[0]],
            "client entries carry timestamps, not ages")

    listing = get("/api/v1/workunits?status=ASSIGNED&limit=100").get_json()
    f.equal(len(listing["workunits"]), 10,
            "filtering by status returns the right count")
    f.check(all(w["status_name"] == "ASSIGNED"
                for w in listing["workunits"]),
            "every returned workunit matches the filter")

    by_client = get("/api/v1/workunits?assigned_to=grvingt-03"
                    "&status=ASSIGNED").get_json()
    f.equal(len(by_client["workunits"]), 4,
            "filtering by client returns the right count")

    by_task = get("/api/v1/workunits?task=polyselect"
                  "&status=ASSIGNED").get_json()
    f.equal([w["wuid"] for w in by_task["workunits"]],
            [EXPECT["stale_from_gone"]],
            "filtering by task returns the right workunit")

    one = get("/api/v1/workunits/" + EXPECT["stale_from_gone"]).get_json()
    f.equal(one["task"], "polyselect", "the task is parsed out of the wuid")
    f.equal(one["identifier"], "5000-5100",
            "the identifier is parsed out of the wuid")
    f.check("workunit" in one, "the workunit body is included")
    f.equal(client.get("/api/v1/workunits/nonexistent",
                       headers=headers).status_code, 404,
            "an unknown workunit is 404")

    log = get("/api/v1/log?tail=5").get_json()
    f.equal(len(log["lines"]), 5, "the log tail honours its bound")
    f.check(log["lines"][-1].endswith("something went wrong"),
            "the log tail returns the *last* lines")

    stats = get("/api/v1/stats").get_json()
    named = {t["name"]: t for t in stats["tasks"]}
    f.check("cputime_las" in named["sieving"].get("times", {}),
            "per-program cpu times are reported")
    f.check(named["polyselect"].get("stats"),
            "the statistics a finished task reported are kept")

    # Polling must be cheap.
    for path in ("/api/v1/progress", "/api/v1/clients",
                 "/api/v1/workunits/summary"):
        first = client.get(path, headers=headers)
        again = client.get(path, headers=dict(
            headers, **{"If-None-Match": first.headers["ETag"]}))
        f.equal(again.status_code, 304, "%s revalidates as 304" % path)
        f.check(first.headers.get("X-Cado-Server-Time"),
                "%s reports the server clock in a header" % path)


# ------------------------------------------------------------------
# actions
# ------------------------------------------------------------------


def section_actions(app, f):
    client = app.test_client()
    headers = authorized(app)

    def counts():
        return client.get("/api/v1/workunits/summary",
                          headers=headers).get_json()["counts"]

    before = counts()

    response = client.post("/api/v1/clients/grvingt-03/reclaim",
                           headers=headers)
    f.equal(response.status_code, 200, "reclaiming a client answers 200")
    result = response.get_json()
    f.equal(len(result["marked"]), EXPECT["reclaimable_from_gone"],
            "the client's current-task workunits are marked")
    f.check(all(w.startswith(NAME + "_sieving_")
                for w in result["marked"]),
            "only sieving workunits were marked")
    f.equal([s["wuid"] for s in result["skipped"]],
            [EXPECT["stale_from_gone"]],
            "the leftover from a finished task is skipped")
    f.check("polyselect" in result["skipped"][0]["reason"],
            "the skip says which task the workunit belongs to")

    after = counts()
    f.equal(after["NEED_RESUBMIT"],
            before["NEED_RESUBMIT"] + EXPECT["reclaimable_from_gone"],
            "the marked workunits are in NEED_RESUBMIT")
    f.equal(after["ASSIGNED"],
            before["ASSIGNED"] - EXPECT["reclaimable_from_gone"],
            "and no longer assigned")
    f.equal(after["CANCELLED"], before["CANCELLED"],
            "nothing was cancelled behind the running task's back")

    # Marking a workunit that belongs to a task which is not running
    # would make ClientServerTask.resubmit_timed_out_wus() charge it to
    # the wrong task, so it must be refused outright.
    response = client.post("/api/v1/workunits/%s/resubmit"
                           % EXPECT["stale_from_gone"], headers=headers)
    f.equal(response.status_code, 409,
            "resubmitting another task's workunit is refused with 409")

    f.equal(client.post("/api/v1/workunits/nonexistent/resubmit",
                        headers=headers).status_code, 404,
            "resubmitting an unknown workunit is 404")

    listing = client.get("/api/v1/workunits?status=ASSIGNED"
                         "&assigned_to=grvingt-01", headers=headers)
    victim = listing.get_json()["workunits"][0]["wuid"]
    f.equal(client.post("/api/v1/workunits/%s/resubmit" % victim,
                        headers=headers).status_code, 200,
            "resubmitting a running task's workunit works")

    # Bulk reclaim, by age.
    result = client.post("/api/v1/workunits/reclaim", headers=headers,
                         json={"older_than": 10 * WUTIMEOUT}).get_json()
    f.equal(result["marked"], [],
            "nothing is old enough for a very long cutoff")

    # Actions need the token too.
    f.equal(client.post("/api/v1/clients/grvingt-01/reclaim").status_code,
            401, "actions are refused without a token")

    # Stop and resume serving, and check a client sees it.
    f.equal(client.post("/api/v1/serving", headers=headers,
                        json={"serving": False}).get_json()["serving"],
            False, "serving can be turned off")
    f.equal(client.get("/workunit", data={"clientid": "t"}).status_code,
            410, "a client is then told to finish")
    f.equal(client.post("/api/v1/serving", headers=headers,
                        json={"serving": True}).get_json()["serving"],
            True, "serving can be turned back on")
    f.equal(client.get("/workunit", data={"clientid": "t"}).status_code,
            200, "and a client gets work again")
    f.equal(client.post("/api/v1/serving", headers=headers,
                        json={}).status_code, 400,
            "a serving request without the field is refused")


# ------------------------------------------------------------------
# doctests
# ------------------------------------------------------------------


def section_doctests(app, f):
    """
    The doctests of the api modules.

    They are run from here rather than added to PYTHON_DOCTEST_SOURCES
    because two of these modules import flask, and this test already
    knows how to skip when flask is absent.
    """
    import doctest
    from cadofactor.api import admin, auth, spec, views

    for module in (spec, auth, views, admin):
        result = doctest.testmod(module)
        f.check(result.failed == 0,
                "%s: %d doctest(s)" % (module.__name__, result.attempted),
                None if not result.failed
                else "%d of %d failed" % (result.failed, result.attempted))
        f.check(result.attempted > 0,
                "%s has doctests at all" % module.__name__)


# ------------------------------------------------------------------
# serve, for test_monitor_cli.sh
# ------------------------------------------------------------------


def section_serve(app, workdir):
    app.serve()
    with open(os.path.join(workdir, "URL"), "w") as f:
        f.write(app.url + "\n")
    print(app.url, flush=True)
    threading.Event().wait()


SECTIONS = {
    "openapi": section_openapi,
    "auth": section_auth,
    "views": section_views,
    "actions": section_actions,
    "doctests": section_doctests,
}


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    section, workdir = argv[1], argv[2]

    import logging
    logging.getLogger().setLevel(logging.CRITICAL)

    app = api_fixture.build(workdir)

    if section == "serve":
        section_serve(app, workdir)
        return 0

    if section not in SECTIONS:
        print("unknown section %r; expected one of %s"
              % (section, ", ".join(sorted(SECTIONS) + ["serve"])))
        return 2

    print("=== %s ===" % section)
    f = Failures()
    SECTIONS[section](app, f)
    print()
    if f.count:
        print("%d check(s) FAILED" % f.count)
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
