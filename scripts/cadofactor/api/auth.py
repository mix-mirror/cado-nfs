"""
Bearer token authentication for the cado-nfs api.

The endpoints that cado-nfs-client.py uses stay unauthenticated. That is
the existing design and it cannot readily change: a client is handed a
url and a certificate fingerprint, never a secret. Those endpoints
remain gated by server.whitelist alone.

Everything else -- the monitoring and administration endpoints under
/api/v1, and the web ui that consumes them -- requires a token which the
server keeps in the working directory, readable only by its owner.
"""

import functools
import hmac
import logging
import os
import secrets
import stat
import time

import flask

# 32 bytes of entropy, rendered by token_urlsafe as 43 characters.
TOKEN_BYTES = 32

# Deliberate cost on a wrong token. The whitelist is the first line of
# defence; this makes guessing pointless even behind a permissive one.
FAILURE_DELAY = 0.25

logger = logging.getLogger("API server")


class TokenError(Exception):
    """
    Raised when an existing token file cannot be trusted.
    """
    pass


def _check_private(path):
    """
    Verify that an existing token file is ours alone.

    Ownership and permission bits are meaningless on non-POSIX
    platforms, where we can only check that the path is a regular file.
    """
    st = os.stat(path)
    if not stat.S_ISREG(st.st_mode):
        raise TokenError("%s is not a regular file" % path)
    if os.name != "posix":
        return
    if st.st_uid != os.getuid():
        raise TokenError("%s is owned by uid %d, not by us (uid %d)"
                         % (path, st.st_uid, os.getuid()))
    if st.st_mode & 0o077:
        raise TokenError("%s is accessible to other users (mode %04o);"
                         " remove it and let the server create a new one"
                         % (path, st.st_mode & 0o7777))


def load_or_create_token(path):
    """
    Return the api token stored at path, creating it if need be.

    An existing token is reused, so that a monitor started against an
    earlier run of the same working directory keeps working. It is
    reused only if it is a regular file owned by us with no access for
    anybody else -- a token that other users can read is not a token, so
    we refuse it loudly rather than pretending it protects anything.

    >>> import tempfile, os
    >>> d = tempfile.mkdtemp()
    >>> p = os.path.join(d, 'c60.api-token')
    >>> t = load_or_create_token(p)
    >>> len(t) >= 40
    True
    >>> oct(os.stat(p).st_mode & 0o777)
    '0o600'
    >>> load_or_create_token(p) == t
    True

    A world-readable token file is rejected rather than trusted:

    >>> os.chmod(p, 0o644)
    >>> try:
    ...     load_or_create_token(p)
    ... except TokenError as e:
    ...     print('accessible to other users' in str(e))
    True

    >>> os.chmod(p, 0o600)
    >>> os.unlink(p)
    >>> os.rmdir(d)
    """
    if os.path.exists(path):
        _check_private(path)
        with open(path, "r") as f:
            token = f.read().strip()
        if token:
            return token
        # An empty file is a leftover from an interrupted create.
        os.unlink(path)

    token = secrets.token_urlsafe(TOKEN_BYTES)
    # O_EXCL so that we never widen the permissions of a file somebody
    # else created in between, and 0600 from the start so that the token
    # is never briefly world-readable.
    fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    try:
        os.write(fd, (token + "\n").encode("ascii"))
    finally:
        os.close(fd)
    return token


def token_from_request(request):
    """
    Extract a bearer token from a request, or None.

    Both spellings are accepted; the header form is what the web ui and
    cado-nfs-monitor.py use.

    >>> class R:
    ...     def __init__(self, h):
    ...         self.headers = h
    >>> token_from_request(R({'Authorization': 'Bearer abc'}))
    'abc'
    >>> token_from_request(R({'Authorization': 'bearer abc'}))
    'abc'
    >>> token_from_request(R({'X-Cado-Token': 'xyz'}))
    'xyz'
    >>> token_from_request(R({'Authorization': 'Basic abc'})) is None
    True
    >>> token_from_request(R({})) is None
    True
    """
    header = request.headers.get("Authorization", "")
    if header[:7].lower() == "bearer ":
        return header[7:].strip()
    direct = request.headers.get("X-Cado-Token")
    if direct:
        return direct.strip()
    return None


def check_token(expected, presented):
    """
    Compare tokens without leaking their contents through timing.

    >>> check_token('secret', 'secret')
    True
    >>> check_token('secret', 'wrong')
    False
    >>> check_token('secret', None)
    False
    >>> check_token(None, 'anything')
    False
    """
    if not expected or not presented:
        return False
    return hmac.compare_digest(str(expected), str(presented))


def require_token(method):
    """
    Decorator gating a flask view on the api token.

    The token is read from the application object, so that the decorator
    can be applied to plain functions as well as to methods of the
    server itself.
    """
    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        expected = getattr(flask.current_app, "api_token", None)
        if expected is None:
            # No token at all means the api is unusable rather than
            # wide open. Say so plainly.
            flask.abort(503, "api token is not configured on the server")
        if not check_token(expected, token_from_request(flask.request)):
            time.sleep(FAILURE_DELAY)
            logger.warning("rejected unauthenticated request for %s"
                           " from %s",
                           flask.request.path, flask.request.remote_addr)
            response = flask.json.jsonify({
                "code": 401,
                "name": "Unauthorized",
                "description":
                    "A valid api token is required. The server wrote one"
                    " to <workdir>/<name>.api-token; send it as"
                    " 'Authorization: Bearer <token>'.",
            })
            response.status_code = 401
            response.headers["WWW-Authenticate"] = 'Bearer realm="cado-nfs"'
            return response
        return method(*args, **kwargs)
    return wrapper


if __name__ == "__main__":
    import doctest
    doctest.testmod()
