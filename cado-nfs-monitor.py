#!/usr/bin/env python3

"""
Watch, and nudge, a running cado-nfs computation from a terminal.

This is a client of the cado-nfs api and nothing else: it never opens
the database, and every single thing it does goes through an endpoint
that is described in /api/v1/openapi.json. That is deliberate. It means
this file can be copied to a laptop on its own, and it means the api
cannot quietly rot while the ui keeps working.

It has no mandatory dependency beyond the standard library. If the rich
package (MIT) happens to be importable, "watch" draws a full-screen
dashboard with it; otherwise the same information is repainted as plain
text. Everything else is unaffected either way.

Typical use, from the machine running cado-nfs.py:

    ./cado-nfs-monitor.py --workdir=/tmp/c120 --server=https://host:8001 \\
            --certsha1=... watch

--workdir is a convenience: it is where the server put both the api
token and its certificate, so giving it saves passing --token and
--certsha1 by hand.
"""

import argparse
import glob
import hashlib
import json
import os
import re
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

try:
    import rich
    import rich.box
    import rich.console
    import rich.layout
    import rich.live
    import rich.panel
    import rich.progress_bar
    import rich.table
    import rich.text
except ImportError:
    rich = None

DEFAULT_TIMEOUT = 20

# Colours, both for rich and for the plain renderer, keyed by the client
# states that /api/v1/clients reports.
STATE_STYLE = {
    "working": ("green", "*"),
    "idle": ("cyan", "-"),
    "stale": ("yellow", "!"),
    "gone": ("red", "X"),
    "unknown": ("dim", "?"),
}

PHASE_MARK = {"done": "+", "running": ">", "pending": " "}


class MonitorError(Exception):
    """
    Anything that should reach the user as a message rather than a
    traceback.
    """
    pass


# --------------------------------------------------------------------
# talking to the server
# --------------------------------------------------------------------


def certificate_fingerprint(pem):
    """
    SHA1 fingerprint of a PEM certificate, as cado-nfs prints it.

    >>> certificate_fingerprint('-----BEGIN CERTIFICATE-----\\n'
    ...                         'AAAA\\n-----END CERTIFICATE-----\\n')
    '29e2dcfbb16f63bb0254df7585a15bb6fb5e927d'
    """
    der = ssl.PEM_cert_to_DER_cert(pem)
    return hashlib.sha1(der).hexdigest()


class Server(object):
    """
    The api server, reached over http or over pinned https.

    We pin the exact certificate by its SHA1 -- which is the only thing
    cado-nfs tells you about it -- and therefore switch hostname
    checking off. The certificate is self-signed and generated for
    whichever address the server happened to bind, so its name is not
    something we can meaningfully verify; the fingerprint is, and it is
    strictly stronger.
    """

    def __init__(self, url, certsha1=None, token=None,
                 timeout=DEFAULT_TIMEOUT):
        self.url = url.rstrip("/")
        self.token = token
        self.timeout = timeout
        self.server_time = None
        self.local_time = None
        parts = urllib.parse.urlparse(self.url)
        self.scheme = parts.scheme
        self.host = parts.hostname
        self.port = parts.port or (443 if self.scheme == "https" else 80)
        self.context = None
        if self.scheme == "https":
            self.context = self._pinned_context(certsha1)

    def _pinned_context(self, certsha1):
        try:
            pem = ssl.get_server_certificate((self.host, self.port))
        except OSError as e:
            raise MonitorError("cannot reach %s: %s" % (self.url, e))
        got = certificate_fingerprint(pem)
        if certsha1 and got.lower() != certsha1.lower().replace(":", ""):
            raise MonitorError(
                "certificate fingerprint mismatch: server presents %s,"
                " --certsha1 says %s. Refusing to continue."
                % (got, certsha1))
        if not certsha1:
            raise MonitorError(
                "this is an https server, so --certsha1 is required."
                " cado-nfs.py prints it when it starts, and it is the"
                " sha1 of %s. The server is presenting %s."
                % ("<workdir>/<name>.server.cert", got))
        context = ssl.create_default_context(cadata=pem)
        # See the class docstring: the fingerprint is the identity.
        context.check_hostname = False
        return context

    def request(self, path, method="GET", body=None):
        url = self.url + path
        data = None
        headers = {"Accept": "application/json"}
        if body is not None:
            data = json.dumps(body).encode("utf-8")
            headers["Content-Type"] = "application/json"
        if self.token:
            headers["Authorization"] = "Bearer " + self.token
        request = urllib.request.Request(url, data=data, headers=headers,
                                         method=method)
        try:
            with urllib.request.urlopen(request, timeout=self.timeout,
                                        context=self.context) as response:
                self._note_server_time(response)
                payload = response.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            self._note_server_time(e)
            raise MonitorError(self._explain(e))
        except urllib.error.URLError as e:
            raise MonitorError("cannot reach %s: %s" % (url, e.reason))
        if not payload:
            return None
        try:
            return json.loads(payload)
        except ValueError:
            raise MonitorError("server did not answer JSON for %s" % path)

    def _note_server_time(self, response):
        """
        Remember the server's clock, which travels in a header rather
        than in the body so that unchanged answers can revalidate.
        Durations are computed against it, never against ours.
        """
        stamp = response.headers.get("X-Cado-Server-Time")
        if stamp:
            try:
                self.server_time = float(stamp)
                self.local_time = time.time()
            except ValueError:
                pass

    def now(self):
        """
        The server's clock, advanced by however long ago we heard it.
        """
        if self.server_time is None:
            return time.time()
        return self.server_time + (time.time() - self.local_time)

    @staticmethod
    def _explain(e):
        detail = ""
        try:
            detail = json.loads(e.read().decode("utf-8")).get(
                "description", "")
        except Exception:
            pass
        if e.code == 401:
            return ("the server rejected our api token."
                    " Pass --workdir, --token-file or --token."
                    + (" (%s)" % detail if detail else ""))
        if e.code == 403:
            return ("the server refused our address."
                    " The monitoring api answers server.ui_whitelist"
                    " only, which is loopback by default; tunnel with"
                    " 'ssh -L' or widen that parameter.")
        return "server returned %d %s%s" % (e.code, e.reason,
                                            ": " + detail if detail else "")


# --------------------------------------------------------------------
# locating the token
# --------------------------------------------------------------------


def find_token(args):
    """
    Work out the api token, in decreasing order of explicitness.

    Passing it as --token is supported but is the worst option: command
    lines are visible to every user on the machine through ps.
    """
    if args.token:
        return args.token
    if args.token_file:
        return read_token_file(args.token_file)
    if os.environ.get("CADO_NFS_API_TOKEN"):
        return os.environ["CADO_NFS_API_TOKEN"]
    if args.workdir:
        candidates = sorted(glob.glob(os.path.join(args.workdir,
                                                   "*.api-token")))
        if len(candidates) == 1:
            return read_token_file(candidates[0])
        if not candidates:
            raise MonitorError("no *.api-token file in %s" % args.workdir)
        raise MonitorError("several token files in %s (%s); say which one"
                           " with --token-file"
                           % (args.workdir,
                              ", ".join(os.path.basename(c)
                                        for c in candidates)))
    return None


def read_token_file(path):
    try:
        with open(path) as f:
            return f.read().strip()
    except OSError as e:
        raise MonitorError("cannot read token file %s: %s" % (path, e))


def find_certsha1(args):
    """
    Derive --certsha1 from the server certificate in the working
    directory, when we have been pointed at one.
    """
    if args.certsha1:
        return None if args.certsha1.lower() == "none" else args.certsha1
    if not args.workdir:
        return None
    candidates = sorted(glob.glob(os.path.join(args.workdir,
                                               "*.server.cert")))
    if len(candidates) != 1:
        return None
    try:
        with open(candidates[0]) as f:
            return certificate_fingerprint(f.read())
    except (OSError, ValueError):
        return None


# --------------------------------------------------------------------
# formatting
# --------------------------------------------------------------------


def human_duration(seconds):
    """
    A duration, at two significant units.

    Short durations keep a decimal: workunits on a small computation
    come back in a fraction of a second, and rounding those to "0s"
    hides exactly the number the reader wanted.

    >>> human_duration(0)
    '0s'
    >>> human_duration(0.42)
    '0.4s'
    >>> human_duration(3.25)
    '3.2s'
    >>> human_duration(45)
    '45s'
    >>> human_duration(3601)
    '1h 0m'
    >>> human_duration(90061)
    '1d 1h'
    >>> human_duration(None)
    '-'
    """
    if seconds is None:
        return "-"
    seconds = max(0, seconds)
    if seconds < 10:
        return ("%.1f" % seconds).rstrip("0").rstrip(".") + "s"
    seconds = int(seconds)
    if seconds < 60:
        return "%ds" % seconds
    if seconds < 3600:
        return "%dm %ds" % (seconds // 60, seconds % 60)
    if seconds < 86400:
        return "%dh %dm" % (seconds // 3600, (seconds % 3600) // 60)
    return "%dd %dh" % (seconds // 86400, (seconds % 86400) // 3600)


def human_number(value):
    """
    A count, abbreviated once it stops being readable.

    >>> human_number(999)
    '999'
    >>> human_number(12345)
    '12.3k'
    >>> human_number(1943221)
    '1.94M'
    >>> human_number(3110000000)
    '3.11G'
    >>> human_number(None)
    '-'
    """
    if value is None:
        return "-"
    value = float(value)
    for limit, suffix in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if abs(value) >= limit:
            return "%.3g%s" % (value / limit, suffix)
    return "%d" % value


def parse_duration(text):
    """
    Parse a human duration such as 30m or 2h into seconds.

    >>> parse_duration('90')
    90.0
    >>> parse_duration('30m')
    1800.0
    >>> parse_duration('2h')
    7200.0
    >>> parse_duration('1d')
    86400.0
    """
    match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*([smhd]?)\s*", text or "")
    if not match:
        raise MonitorError("cannot read %r as a duration; try 30m, 2h, 1d"
                           % text)
    scale = {"": 1, "s": 1, "m": 60, "h": 3600, "d": 86400}
    return float(match.group(1)) * scale[match.group(2)]


def bar(fraction, width=24):
    """
    An ASCII progress bar.

    >>> bar(0.5, 10)
    '[#####-----]'
    >>> bar(0, 4)
    '[----]'
    >>> bar(1, 4)
    '[####]'
    >>> bar(None, 4)
    '[????]'
    """
    if fraction is None:
        return "[" + "?" * width + "]"
    filled = int(round(max(0.0, min(1.0, fraction)) * width))
    return "[" + "#" * filled + "-" * (width - filled) + "]"


def task_line(task, width=24):
    """
    One line of the pipeline, as the plain renderer prints it.
    """
    mark = PHASE_MARK.get(task["phase"], " ")
    title = task.get("title") or task["name"]
    if task["phase"] == "running" and task.get("achievement") is not None:
        achievement = task["achievement"]
        return " %s %-46s %6.2f%% %s" % (mark, title[:46],
                                         100.0 * achievement,
                                         bar(achievement, width))
    if task["phase"] == "done":
        return " %s %-46s %s" % (mark, title[:46], "done")
    return " %s %s" % (mark, title[:46])


def highlight_line(task):
    """
    The task-specific numbers worth showing under the bar.
    """
    h = task.get("highlights") or {}
    bits = []
    if "rels_found" in h and "rels_wanted" in h:
        bits.append("relations %s / %s" % (human_number(h["rels_found"]),
                                           human_number(h["rels_wanted"])))
    if "qnext" in h:
        bits.append("q at %s" % human_number(h["qnext"]))
    if "adnext" in h:
        bits.append("ad at %s" % human_number(h["adnext"]))
    if task.get("eta"):
        bits.append("ETA %s" % task["eta"])
    return "   " + ",  ".join(bits) if bits else ""


# --------------------------------------------------------------------
# gathering
# --------------------------------------------------------------------


def gather(server, want_log=0):
    """
    One round of everything the dashboard shows.
    """
    snapshot = {
        "info": server.request("/api/v1/info"),
        "progress": server.request("/api/v1/progress"),
        "summary": server.request("/api/v1/workunits/summary"),
        "clients": server.request("/api/v1/clients"),
        "now": None,
    }
    if want_log:
        snapshot["log"] = server.request("/api/v1/log?tail=%d" % want_log)
    snapshot["now"] = server.now()
    return snapshot


def current_task(progress):
    for task in progress.get("tasks", []):
        if task["phase"] == "running":
            return task
    return None


# --------------------------------------------------------------------
# plain rendering
# --------------------------------------------------------------------


def render_plain(snapshot, verbose=True):
    info = snapshot["info"]
    progress = snapshot["progress"]
    summary = snapshot["summary"]
    clients = snapshot["clients"]
    now = snapshot["now"]
    out = []

    kind = (info.get("computation_desc")
            or info.get("computation") or "computation")
    header = "cado-nfs %s -- %s" % (kind, info.get("name") or "?")
    if info.get("starttime"):
        header += "   running for %s" % human_duration(
            now - float(info["starttime"]))
    if progress.get("finished"):
        header += "   [FINISHED]"
    elif not info.get("serving_workunits", True):
        header += "   [NOT SERVING WORKUNITS]"
    out.append(header)
    out.append("=" * max(len(header), 60))
    out.append("")

    for task in progress.get("tasks", []):
        if not verbose and task["phase"] == "pending":
            continue
        out.append(task_line(task))
        if task["phase"] == "running":
            line = highlight_line(task)
            if line.strip():
                out.append(line)
    out.append("")

    counts = summary.get("counts", {})
    out.append("Workunits   available %-6d assigned %-6d resubmit %-6d"
               % (counts.get("AVAILABLE", 0), counts.get("ASSIGNED", 0),
                  counts.get("NEED_RESUBMIT", 0)))
    out.append("            done %-11d failed %-8d total %d"
               % (counts.get("VERIFIED_OK", 0)
                  + counts.get("RECEIVED_OK", 0),
                  counts.get("VERIFIED_ERROR", 0)
                  + counts.get("RECEIVED_ERROR", 0),
                  summary.get("total", 0)))
    out.append("")

    rows = clients.get("clients", [])
    tally = ", ".join("%d %s" % (n, s)
                      for s, n in sorted(clients.get("counts", {}).items()))
    out.append("Clients (%d): %s" % (len(rows), tally or "none yet"))
    if rows:
        out.append("  %-24s %-8s %7s %7s %7s %9s"
                   % ("client", "state", "flight", "done", "failed",
                      "last seen"))
        for client in rows:
            last = client.get("last_seen")
            out.append("  %s %-22s %-8s %7d %7d %7d %9s"
                       % (STATE_STYLE.get(client["state"], ("", "?"))[1],
                          client["clientid"][:22],
                          client["state"],
                          client["in_flight"], client["completed"],
                          client["failed"],
                          human_duration(None if last is None
                                         else now - last)))
    stale = [c for c in rows if c["state"] in ("stale", "gone")
             and c["in_flight"]]
    if stale:
        out.append("")
        out.append("  %d client(s) hold %d workunit(s) but have gone"
                   " quiet. Reclaim with:" %
                   (len(stale), sum(c["in_flight"] for c in stale)))
        for client in stale:
            out.append("      cado-nfs-monitor.py ... clients reclaim %s"
                       % client["clientid"])

    if "log" in snapshot:
        out.append("")
        out.append("Log:")
        for line in snapshot["log"].get("lines", []):
            out.append("  " + line)
    return "\n".join(out)


# --------------------------------------------------------------------
# rich rendering
# --------------------------------------------------------------------


def render_rich(snapshot):
    """
    Build the renderable that "watch" displays when rich is available.
    """
    info = snapshot["info"]
    progress = snapshot["progress"]
    summary = snapshot["summary"]
    clients = snapshot["clients"]
    now = snapshot["now"]

    title = rich.text.Text()
    title.append("cado-nfs ", style="bold")
    title.append(str(info.get("computation_desc")
                     or info.get("computation") or "computation"))
    title.append("  %s" % (info.get("name") or "?"), style="bold cyan")
    if info.get("starttime"):
        title.append("   running for %s"
                     % human_duration(now - float(info["starttime"])),
                     style="dim")
    if progress.get("finished"):
        title.append("   FINISHED", style="bold green")
    elif not info.get("serving_workunits", True):
        title.append("   NOT SERVING WORKUNITS", style="bold yellow")

    pipeline = rich.table.Table.grid(padding=(0, 1), expand=True)
    pipeline.add_column(width=2)
    # One row per task and no wrapping, so that the panel is exactly as
    # tall as the layout below reserves for it.
    pipeline.add_column(ratio=1, no_wrap=True, overflow="ellipsis")
    pipeline.add_column(width=26)
    pipeline.add_column(width=8, justify="right")
    for task in progress.get("tasks", []):
        phase = task["phase"]
        if phase == "done":
            mark, style, right = "+", "green", "done"
        elif phase == "running":
            mark, style, right = ">", "bold white", ""
        else:
            mark, style, right = " ", "dim", ""
        achievement = task.get("achievement")
        if phase == "running" and achievement is not None:
            widget = rich.progress_bar.ProgressBar(
                total=1.0, completed=max(0.0, min(1.0, achievement)),
                width=26)
            right = "%.2f%%" % (100.0 * achievement)
        elif phase == "done":
            widget = rich.text.Text("")
        else:
            widget = rich.text.Text("")
        pipeline.add_row(rich.text.Text(mark, style=style),
                         rich.text.Text(task.get("title") or task["name"],
                                        style=style),
                         widget,
                         rich.text.Text(right, style=style))

    running = current_task(progress)
    subtitle = highlight_line(running).strip() if running else ""
    pipeline_panel = rich.panel.Panel(
        pipeline, title="pipeline", subtitle=subtitle or None,
        border_style="blue", box=rich.box.ROUNDED)

    counts = summary.get("counts", {})
    wu = rich.table.Table.grid(padding=(0, 2))
    wu.add_column(justify="right", style="dim")
    wu.add_column(justify="right")
    for label, value, style in (
            ("available", counts.get("AVAILABLE", 0), "cyan"),
            ("assigned", counts.get("ASSIGNED", 0), "green"),
            ("to resubmit", counts.get("NEED_RESUBMIT", 0), "yellow"),
            ("done", counts.get("VERIFIED_OK", 0)
             + counts.get("RECEIVED_OK", 0), "white"),
            ("failed", counts.get("VERIFIED_ERROR", 0)
             + counts.get("RECEIVED_ERROR", 0), "red"),
            ("total", summary.get("total", 0), "bold")):
        wu.add_row(label, rich.text.Text(str(value), style=style))
    wu_panel = rich.panel.Panel(wu, title="workunits",
                                border_style="blue", box=rich.box.ROUNDED)

    table = rich.table.Table(box=rich.box.SIMPLE, expand=True,
                             pad_edge=False)
    table.add_column("client", ratio=1, no_wrap=True)
    table.add_column("state", width=8)
    table.add_column("flight", justify="right", width=6)
    table.add_column("done", justify="right", width=7)
    table.add_column("failed", justify="right", width=6)
    table.add_column("share", justify="right", width=6)
    table.add_column("last seen", justify="right", width=10)
    table.add_column("pace", justify="right", width=9)
    for client in clients.get("clients", []):
        style = STATE_STYLE.get(client["state"], ("dim", "?"))[0]
        last = client.get("last_seen")
        table.add_row(client["clientid"],
                      rich.text.Text(client["state"], style=style),
                      str(client["in_flight"]),
                      str(client["completed"]),
                      rich.text.Text(str(client["failed"]),
                                     style="red" if client["failed"]
                                     else "dim"),
                      "%.0f%%" % (100.0 * client.get("share", 0.0)),
                      human_duration(None if last is None
                                     else now - last),
                      rich.text.Text(
                          human_duration(client.get("typical_turnaround")),
                          style="dim" if client.get("typical_turnaround")
                          is None else ""))
    tally = "  ".join("%d %s" % (n, s) for s, n
                      in sorted(clients.get("counts", {}).items()))
    clients_panel = rich.panel.Panel(
        table, title="clients", subtitle=tally or None,
        border_style="blue", box=rich.box.ROUNDED)

    layout = rich.layout.Layout()
    layout.split_column(
        rich.layout.Layout(rich.panel.Panel(title, box=rich.box.HEAVY,
                                            border_style="cyan"),
                           size=3),
        rich.layout.Layout(name="middle", size=len(
            progress.get("tasks", [])) + 2),
        rich.layout.Layout(name="bottom"))
    layout["middle"].split_row(
        rich.layout.Layout(pipeline_panel, ratio=3),
        rich.layout.Layout(wu_panel, ratio=1))
    layout["bottom"].update(clients_panel)
    return layout


# --------------------------------------------------------------------
# commands
# --------------------------------------------------------------------


def emit(args, payload, text):
    if args.json:
        json.dump(payload, sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
    else:
        print(text)


def cmd_status(server, args):
    snapshot = gather(server)
    emit(args, snapshot, render_plain(snapshot, verbose=not args.brief))
    return 0


def cmd_watch(server, args):
    if args.json:
        raise MonitorError("--json makes no sense with watch;"
                           " use status --json in a loop")
    if rich is not None and sys.stdout.isatty():
        return watch_rich(server, args)
    return watch_plain(server, args)


def watch_rich(server, args):
    console = rich.console.Console()
    with rich.live.Live(console=console, screen=True,
                        refresh_per_second=4) as live:
        while True:
            try:
                live.update(render_rich(gather(server)))
            except MonitorError as e:
                live.update(rich.panel.Panel(
                    rich.text.Text(str(e), style="red"),
                    title="lost contact with the server",
                    border_style="red"))
            try:
                time.sleep(args.interval)
            except KeyboardInterrupt:
                return 0


def watch_plain(server, args):
    while True:
        try:
            text = render_plain(gather(server))
        except MonitorError as e:
            text = "lost contact with the server: %s" % e
        if sys.stdout.isatty():
            # Home and clear-to-end, rather than a full clear, so the
            # display does not flicker.
            sys.stdout.write("\033[H\033[J")
        sys.stdout.write(text + "\n")
        sys.stdout.flush()
        try:
            time.sleep(args.interval)
        except KeyboardInterrupt:
            return 0


def cmd_clients(server, args):
    if args.reclaim:
        result = server.request("/api/v1/clients/%s/reclaim"
                                % urllib.parse.quote(args.reclaim, safe=""),
                                method="POST")
        return report_action(args, result)
    payload = server.request("/api/v1/clients")
    now = server.now()
    lines = ["%-24s %-8s %7s %7s %7s %6s %10s %10s"
             % ("client", "state", "flight", "done", "failed", "share",
                "last seen", "pace")]
    for client in payload.get("clients", []):
        last = client.get("last_seen")
        lines.append("%-24s %-8s %7d %7d %7d %5.0f%% %10s %10s"
                     % (client["clientid"][:24], client["state"],
                        client["in_flight"], client["completed"],
                        client["failed"],
                        100.0 * client.get("share", 0.0),
                        human_duration(None if last is None
                                       else now - last),
                        human_duration(client.get("typical_turnaround"))))
    stale = [c for c in payload.get("clients", [])
             if c["state"] in ("stale", "gone")]
    if stale:
        lines.append("")
        lines.append("'pace' is how long this client's workunits have"
                     " recently been taking, and")
        lines.append("is what its staleness threshold is derived from:")
        for c in stale:
            lines.append("  %-24s silent for %s, counted stale after %s"
                         " (%s)"
                         % (c["clientid"][:24],
                            human_duration(None if c.get("last_seen") is None
                                           else now - c["last_seen"]),
                            human_duration(c.get("stale_after")),
                            "its own pace"
                            if c.get("liveness_basis") == "turnaround"
                            else "tasks.wutimeout; too few samples yet"))
    emit(args, payload, "\n".join(lines))
    return 0


def cmd_wu_list(server, args):
    query = {"limit": args.limit}
    if args.status:
        query["status"] = args.status
    if args.assigned_to:
        query["assigned_to"] = args.assigned_to
    if args.result_from:
        query["result_from"] = args.result_from
    if args.task:
        query["task"] = args.task
    if args.older_than:
        query["assigned_older_than"] = int(parse_duration(args.older_than))
    payload = server.request("/api/v1/workunits?"
                             + urllib.parse.urlencode(query))
    now = server.now()
    lines = ["%-44s %-14s %-20s %10s"
             % ("workunit", "status", "client", "age")]
    for wu in payload.get("workunits", []):
        stamp = wu.get("timeresult") or wu.get("timeassigned") \
            or wu.get("timecreated")
        lines.append("%-44s %-14s %-20s %10s"
                     % (wu["wuid"][:44], wu["status_name"],
                        (wu.get("assignedclient")
                         or wu.get("resultclient") or "-")[:20],
                        human_duration(None if stamp is None
                                       else now - stamp)))
    emit(args, payload, "\n".join(lines))
    return 0


def cmd_wu_show(server, args):
    payload = server.request("/api/v1/workunits/%s"
                             % urllib.parse.quote(args.wuid, safe=""))
    now = server.now()
    lines = []
    for key in ("wuid", "status_name", "task", "identifier", "attempt",
                "assignedclient", "resultclient", "errorcode"):
        if payload.get(key) is not None:
            lines.append("%-16s %s" % (key, payload[key]))
    for key in ("timecreated", "timeassigned", "timeresult",
                "timeverified"):
        if payload.get(key):
            lines.append("%-16s %s ago" % (key,
                                           human_duration(now
                                                          - payload[key])))
    for f in payload.get("files") or []:
        lines.append("%-16s %s" % ("file", f.get("path")))
    emit(args, payload, "\n".join(lines))
    return 0


def cmd_wu_resubmit(server, args):
    result = server.request("/api/v1/workunits/%s/resubmit"
                            % urllib.parse.quote(args.wuid, safe=""),
                            method="POST")
    return report_action(args, result)


def cmd_wu_reclaim(server, args):
    body = {}
    if args.older_than:
        body["older_than"] = int(parse_duration(args.older_than))
    result = server.request("/api/v1/workunits/reclaim", method="POST",
                            body=body)
    return report_action(args, result)


def cmd_log(server, args):
    seen = 0
    while True:
        payload = server.request("/api/v1/log?tail=%d" % args.tail)
        lines = payload.get("lines", [])
        if args.json:
            emit(args, payload, "")
            return 0
        for line in lines[seen:] if seen else lines:
            print(line)
        if not args.follow:
            return 0
        seen = len(lines)
        try:
            time.sleep(args.interval)
        except KeyboardInterrupt:
            return 0


def cmd_serving(server, args):
    if args.state is None:
        payload = server.request("/api/v1/serving")
    else:
        payload = server.request("/api/v1/serving", method="POST",
                                 body={"serving": args.state == "on"})
    emit(args, payload,
         "serving workunits: %s" % ("yes" if payload["serving"] else "no"))
    return 0


def cmd_raw(server, args):
    payload = server.request(args.path, method=args.method)
    json.dump(payload, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


def report_action(args, result):
    if args.json:
        json.dump(result, sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
        return 0
    print(result.get("message", ""))
    for wuid in result.get("marked", []):
        print("  marked  %s" % wuid)
    for entry in result.get("skipped", []):
        print("  skipped %s -- %s" % (entry.get("wuid"),
                                      entry.get("reason")))
    return 0 if result.get("marked") else 1


# --------------------------------------------------------------------
# command line
# --------------------------------------------------------------------


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__.strip().split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="The api this drives is described at"
               " <server>/api/v1/openapi.json, and rendered at"
               " <server>/api/docs.")
    parser.add_argument("--server", required=True,
                        help="url of the cado-nfs api server, as"
                             " cado-nfs.py prints it")
    parser.add_argument("--certsha1",
                        help="sha1 of the server certificate, as"
                             " cado-nfs.py prints it; 'None' for a plain"
                             " http server. Deduced from --workdir when"
                             " that is given.")
    parser.add_argument("--workdir",
                        help="the computation's working directory, used"
                             " to find the api token and the server"
                             " certificate")
    parser.add_argument("--token-file",
                        help="file holding the api token")
    parser.add_argument("--token",
                        help="the api token itself. Prefer --workdir or"
                             " --token-file: a command line is visible"
                             " to every user on the machine.")
    parser.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT,
                        help="seconds to wait for the server"
                             " (default: %(default)s)")
    parser.add_argument("--json", action="store_true",
                        help="print the api's answer verbatim, for"
                             " scripting")

    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("status", help="one-shot summary")
    p.add_argument("--brief", action="store_true",
                   help="hide tasks that have not started")
    p.set_defaults(func=cmd_status)

    p = sub.add_parser("watch", help="live dashboard")
    p.add_argument("--interval", type=float, default=3.0,
                   help="seconds between refreshes (default: %(default)s)")
    p.set_defaults(func=cmd_watch)

    p = sub.add_parser("clients", help="what each client is up to")
    p.add_argument("--reclaim", metavar="CLIENTID",
                   help="put the workunits this client holds back in the"
                        " pool")
    p.set_defaults(func=cmd_clients)

    p = sub.add_parser("wu", help="inspect and requeue workunits")
    wusub = p.add_subparsers(dest="wucommand", required=True)

    q = wusub.add_parser("list", help="browse the workunits table")
    q.add_argument("--status", help="status name or number")
    q.add_argument("--assigned-to", metavar="CLIENTID")
    q.add_argument("--result-from", metavar="CLIENTID")
    q.add_argument("--task")
    q.add_argument("--older-than", metavar="DURATION",
                   help="only workunits assigned longer ago than this,"
                        " e.g. 30m")
    q.add_argument("--limit", type=int, default=40)
    q.set_defaults(func=cmd_wu_list)

    q = wusub.add_parser("show", help="one workunit in detail")
    q.add_argument("wuid")
    q.set_defaults(func=cmd_wu_show)

    q = wusub.add_parser("resubmit", help="put one workunit back in the"
                                          " pool")
    q.add_argument("wuid")
    q.set_defaults(func=cmd_wu_resubmit)

    q = wusub.add_parser("reclaim",
                         help="put every long-outstanding workunit back"
                              " in the pool")
    q.add_argument("--older-than", metavar="DURATION",
                   help="age of the assignment; defaults to the"
                        " computation's tasks.wutimeout")
    q.set_defaults(func=cmd_wu_reclaim)

    p = sub.add_parser("log", help="tail the computation's log")
    p.add_argument("--tail", type=int, default=40)
    p.add_argument("--follow", "-f", action="store_true")
    p.add_argument("--interval", type=float, default=3.0)
    p.set_defaults(func=cmd_log)

    p = sub.add_parser("serving",
                       help="whether the server hands workunits out")
    p.add_argument("state", nargs="?", choices=["on", "off"])
    p.set_defaults(func=cmd_serving)

    p = sub.add_parser("raw", help="call an endpoint and print the JSON")
    p.add_argument("path", help="e.g. /api/v1/stats")
    p.add_argument("--method", default="GET")
    p.set_defaults(func=cmd_raw)

    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        server = Server(args.server,
                        certsha1=find_certsha1(args),
                        token=find_token(args),
                        timeout=args.timeout)
        return args.func(server, args)
    except MonitorError as e:
        print("%s: %s" % (os.path.basename(sys.argv[0]), e),
              file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        return 0
    except BrokenPipeError:
        # e.g. piped into head(1)
        try:
            sys.stdout.close()
        except Exception:
            pass
        return 0


if __name__ == "__main__":
    if "--doctest" in sys.argv:
        import doctest
        sys.exit(bool(doctest.testmod().failed))
    sys.exit(main())
