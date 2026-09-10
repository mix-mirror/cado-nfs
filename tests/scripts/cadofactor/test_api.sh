#!/usr/bin/env bash

# Run one section of the monitoring api checks. See test_api.py.

set -e

: ${wdir:?missing}
: ${CADO_NFS_SOURCE_DIR:?missing}

section="$1"

if [ -z "$section" ] ; then
    echo "usage: $0 <section>" >&2
    exit 1
fi

export PYTHONPATH="${CADO_NFS_SOURCE_DIR}/scripts${PYTHONPATH:+:$PYTHONPATH}"

# The api server needs flask. It is an optional dependency of cado-nfs
# as a whole, so skip rather than fail when it is missing -- exactly as
# a build without it would behave.
if ! python3 -c 'import flask' 2>/dev/null ; then
    echo "flask is not installed, skipping the api tests"
    exit 0
fi

exec python3 "${CADO_NFS_SOURCE_DIR}/tests/scripts/cadofactor/test_api.py" \
    "$section" "$wdir/api-$section"
