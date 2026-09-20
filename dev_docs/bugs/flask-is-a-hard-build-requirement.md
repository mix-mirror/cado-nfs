# flask is required to configure the build, not merely to run the server

## Summary

`config/python.cmake` lists `flask` in `PYTHON_REQUIRED_MODS`, so cmake
refuses to configure at all when it is absent -- even for a build that
will never start a server. The api server, the dashboard and the
monitoring endpoints are all things one can do without; the build is
not.

## Reproducer

On a machine whose python has no flask (a fresh Grid'5000 node, for
instance):

    $ make
    CMake Error at config/python.cmake:49 (message):
      Importing the flask Python module failed.  This may be caused by the
      flask library package missing on your system. [...]
    -- Configuring incomplete, errors occurred!

The message goes on to suggest `./scripts/setup-venv.sh`, which does
work, but the point is that the build stops.

## Suggested fix

Move `flask` from `PYTHON_REQUIRED_MODS` to `PYTHON_OPTIONAL_MODS`.
`HAVE_PYTHON_FLASK` is already set to 0 in that case, so the
information is available to whatever needs to react to it. The tests
that need a server already skip themselves when flask is missing --
`tests/scripts/cadofactor/test_api.sh` and `test_monitor_cli.sh` both
start with an import check -- so making it optional costs those tests
nothing beyond what they already handle.

`requests` is in the same list and has a stronger claim to being
required, since `cado-nfs-client.py` cannot run without it; but that
too is a runtime need of one script rather than a build-time need.
