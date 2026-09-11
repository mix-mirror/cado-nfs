#!/usr/bin/env bash

# Drive cado-nfs-monitor.py against a real http server, so that the
# whole path -- argument parsing, urllib, token discovery, the api
# itself -- is exercised rather than only the flask test client.

set -e

: ${wdir:?missing}
: ${CADO_NFS_SOURCE_DIR:?missing}

export PYTHONPATH="${CADO_NFS_SOURCE_DIR}/scripts${PYTHONPATH:+:$PYTHONPATH}"

if ! python3 -c 'import flask' 2>/dev/null ; then
    echo "flask is not installed, skipping the monitor tests"
    exit 0
fi

server_wdir="$wdir/monitor"
mkdir -p "$server_wdir"

python3 "${CADO_NFS_SOURCE_DIR}/tests/scripts/cadofactor/test_api.py" \
    serve "$server_wdir" > "$wdir/server.out" 2>"$wdir/server.err" &
server_pid=$!

cleanup() {
    kill "$server_pid" 2>/dev/null || true
    wait "$server_pid" 2>/dev/null || true
}
trap cleanup EXIT

# A busy runner is much slower than a developer's machine: allow well
# over the nominal startup, but in whole seconds, since not every sh
# has a sleep that takes fractions.
for _ in $(seq 75) ; do
    [ -s "$server_wdir/URL" ] && break
    sleep 1
done

if ! [ -s "$server_wdir/URL" ] ; then
    echo "the test server never came up" >&2
    echo "--- its stdout ---" >&2
    cat "$wdir/server.out" >&2 || :
    echo "--- its stderr ---" >&2
    cat "$wdir/server.err" >&2 || :
    exit 1
fi

url=$(cat "$server_wdir/URL")
monitor=("python3" "${CADO_NFS_SOURCE_DIR}/cado-nfs-monitor.py"
         "--server=$url" "--workdir=$server_wdir" "--certsha1=None")

set -x

# Its own doctests, which cover the formatting helpers.
python3 "${CADO_NFS_SOURCE_DIR}/cado-nfs-monitor.py" --doctest

# The token is found in the working directory, and every read-only
# subcommand answers.
"${monitor[@]}" status > "$wdir/status.txt"
grep -q "Lattice Sieving" "$wdir/status.txt"
grep -q "alpha-03" "$wdir/status.txt"

"${monitor[@]}" clients | grep -q "alpha-01"
"${monitor[@]}" wu list --status=ASSIGNED --limit=5 | grep -q ASSIGNED
"${monitor[@]}" wu show c60_sieving_970000-971000 | grep -q alpha-02
"${monitor[@]}" log --tail=5 | grep -q "Lattice Sieving"
"${monitor[@]}" serving | grep -q "serving workunits: yes"

# --json must be machine-readable on every subcommand that has one.
for subcommand in status clients ; do
    "${monitor[@]}" --json "$subcommand" | python3 -m json.tool > /dev/null
done

# An action, and its effect.
"${monitor[@]}" clients --reclaim alpha-03 | grep -q "marked for resubmission"
"${monitor[@]}" wu list --status=NEED_RESUBMIT | grep -q NEED_RESUBMIT

# A workunit from a task that is not running must be refused.
if "${monitor[@]}" wu resubmit c60_polyselect_5000-5100 ; then
    echo "resubmitting another task's workunit should have failed" >&2
    exit 1
fi

# Without a token, the monitor must say so rather than crash.
if python3 "${CADO_NFS_SOURCE_DIR}/cado-nfs-monitor.py" \
        "--server=$url" --certsha1=None status 2>"$wdir/notoken.txt" ; then
    echo "the monitor should have refused to work without a token" >&2
    exit 1
fi
grep -q "token" "$wdir/notoken.txt"

set +x
echo "monitor cli checks passed"
