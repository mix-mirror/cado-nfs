#!/usr/bin/env bash

# Check that every route the api server registers appears in the
# OpenAPI document it serves.
#
# This is a consistency check on the source, not a test of behaviour:
# its answer cannot differ between one machine and another, so it needs
# to run once, not on every platform in CI. Hence a commit hook rather
# than a ctest entry.
#
# Like scripts/check_python_lint.sh, it does nothing unless the commit
# actually touches the files it is about. And like the api itself, it
# treats flask as optional: without it there is nothing to check, and
# saying so is better than failing.

set -e

# A commit hook must leave no trace: without this, importing the api
# writes .pyc files into the source tree, which the very next check in
# check_repo_policies.sh (check_file_lists.pl) then objects to.
export PYTHONDONTWRITEBYTECODE=1

changed_files() {
    if [ "$1" = "--all" ] ; then
        git ls-files 'scripts/cadofactor/api/*.py' \
                     'scripts/cadofactor/api_server.py'
    elif [ "$*" ] ; then
        echo "$@"
    else
        git diff --cached --name-only --diff-filter=ACM \
            'scripts/cadofactor/api/*.py' \
            'scripts/cadofactor/api_server.py'
    fi
}

files=(`changed_files "$@"`)

if [ "${#files[@]}" = 0 ] ; then
    exit 0
fi

here=$(cd "$(dirname "$0")" && pwd)
top=$(dirname "$here")

if ! PYTHONPATH="$top/scripts${PYTHONPATH:+:$PYTHONPATH}" \
        python3 -c 'import flask' 2>/dev/null ; then
    echo "check_openapi_coverage: flask is absent, nothing to check"
    exit 0
fi

PYTHONPATH="$top/scripts${PYTHONPATH:+:$PYTHONPATH}" python3 - "$top" <<'EOF'
import logging
import os
import sys
import tempfile

logging.disable(logging.CRITICAL)

top = sys.argv[1]

from cadofactor import wudb                                # noqa: E402
from cadofactor.api_server import ApiServer                # noqa: E402
from cadofactor.api.spec import rule_to_openapi_path       # noqa: E402

problems = []

with tempfile.TemporaryDirectory() as workdir:
    db = wudb.DBFactory("db:sqlite3://%s/check.db" % workdir, create=True)
    wudb.WuAccess(db.connect()).create_tables()
    app = ApiServer(None, 0, db, threaded=False,
                    uploaddir=os.path.join(workdir, "upload"),
                    nrsubdir=0, cafile=None,
                    whitelist=["127.0.0.1/32"],
                    workdir=workdir, name="check")
    document = app.openapi_spec()

    reachable = {rule_to_openapi_path(str(rule))
                 for rule in app.url_map.iter_rules()
                 if rule.endpoint != "static"}
    documented = set(document["paths"])

    for path in sorted(reachable - documented):
        problems.append("%s is reachable but appears in no"
                        " OpenAPI document" % path)
    for path in sorted(documented - reachable):
        problems.append("%s is documented but is not reachable" % path)

    for path, entry in sorted(document["paths"].items()):
        for method, operation in sorted(entry.items()):
            for field in ("summary", "operationId", "responses", "tags"):
                if not operation.get(field):
                    problems.append("%s %s has no %s"
                                    % (method.upper(), path, field))

    ids = [op["operationId"] for e in document["paths"].values()
           for op in e.values()]
    for name in sorted({i for i in ids if ids.count(i) > 1}):
        problems.append("operationId %s is used more than once" % name)

if problems:
    print("The api and its OpenAPI document disagree:", file=sys.stderr)
    for problem in problems:
        print("  %s" % problem, file=sys.stderr)
    print("", file=sys.stderr)
    print("Routes and documentation come from the same @api_route"
          " declaration, so this normally means a route was registered"
          " by hand. See scripts/cadofactor/api/spec.py.", file=sys.stderr)
    sys.exit(1)
EOF
