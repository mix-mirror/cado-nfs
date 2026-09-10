"""
A pool of database sessions for the api server.

The server used to keep one connection per thread, in a dictionary
keyed by thread identity. That is fine while the server is
single-threaded, which is what server.threaded=False gives you, and it
is quietly disastrous the moment it is not: werkzeug's threaded server
is a ThreadingMixIn, so it starts *a new thread for every request*.
Keyed by thread identity, that means

 - a fresh sqlite connection per request, and
 - a fresh DictDbDirectAccess per request, whose constructor issues
   CREATE TABLE IF NOT EXISTS inside an EXCLUSIVE transaction -- so
   every request would take a write lock on a database that the
   computation itself is trying to use, and

 - a dictionary that grows for as long as the server runs, since dead
   threads are never removed from it.

So sessions are pooled and reused instead. A request borrows one and
gives it back; the expensive setup happens once per pooled session
rather than once per request. The pool is bounded, which also bounds
how many writers can queue up against sqlite at once -- there is
nothing to be gained by letting a thousand threads contend for a lock
that only one of them can hold.

Sessions must not be used by two threads at once. Borrowing guarantees
that, which is what makes it safe to open their connections with
DBFactory.connect(shared_across_threads=True) -- see DbSession.
"""

import logging
import queue
import threading

from cadofactor import wudb
from cadofactor.database import DictDbDirectAccess
from cadofactor.database.base import conn_close
from cadofactor.api.views import (SERVER_STATE_TABLE, API_OVERRIDES_TABLE,
                                  CLIENT_INFO_TABLE)

logger = logging.getLogger("API server")

REGISTERED_FILENAMES_TABLE = "server_registered_filenames"

# How long to wait for a session before deciding that something is
# wrong. Requests queue behind the pool by design; they should not
# queue forever.
BORROW_TIMEOUT = 120.0


class DbSession(object):
    """
    One database connection, and the handful of objects built on it.

    Everything here is created once, when the session is created. Each
    of these constructors takes an EXCLUSIVE transaction to create its
    table if it is missing, which is exactly the cost this class exists
    to stop paying per request.

    On threads: a pooled session outlives the request that created it,
    so it will sooner or later be borrowed by a different thread than
    the one that opened its connection. sqlite3 forbids that by
    default, and rightly. The caller therefore opens the connection
    with DBFactory.connect(shared_across_threads=True), which lifts
    that specific guard -- permissible only because the pool hands a
    session to exactly one borrower at a time, which is the condition
    the guard exists to enforce.
    """

    def __init__(self, connect):
        self.connection = connect()
        self.wuaccess = wudb.WuAccess(self.connection)
        self.registered_filenames = DictDbDirectAccess(
            self.connection, REGISTERED_FILENAMES_TABLE)
        self.server_state = DictDbDirectAccess(self.connection,
                                               SERVER_STATE_TABLE)
        self.overrides = DictDbDirectAccess(self.connection,
                                            API_OVERRIDES_TABLE)
        self.client_info = DictDbDirectAccess(self.connection,
                                              CLIENT_INFO_TABLE)

    def close(self):
        try:
            if callable(conn_close):
                conn_close(self.connection)
            else:
                self.connection.close()
        except Exception as e:
            logger.debug("while closing a db session: %s", e)


class DbSessionPool(object):
    """
    A bounded pool of DbSession objects.

    >>> made = []
    >>> class FakeSession:
    ...     def __init__(self): made.append(self); self.closed = False
    ...     def close(self): self.closed = True
    >>> pool = DbSessionPool(FakeSession, maxsize=2)
    >>> a = pool.borrow()
    >>> b = pool.borrow()
    >>> len(made)
    2
    >>> pool.release(a)
    >>> c = pool.borrow()          # reuses the one just released
    >>> len(made)
    2
    >>> c is a
    True
    >>> pool.release(b); pool.release(c)
    >>> pool.close()
    >>> [s.closed for s in made]
    [True, True]
    """

    def __init__(self, factory, maxsize):
        self._factory = factory
        self._maxsize = max(1, int(maxsize))
        self._idle = queue.LifoQueue()
        self._lock = threading.Lock()
        self._created = 0
        self._closed = False

    def borrow(self):
        """
        Take a session, creating one if the pool is not yet full and
        waiting for one to come back if it is.
        """
        try:
            return self._idle.get_nowait()
        except queue.Empty:
            pass
        with self._lock:
            if self._created < self._maxsize and not self._closed:
                self._created += 1
                create = True
            else:
                create = False
        if create:
            try:
                return self._factory()
            except Exception:
                with self._lock:
                    self._created -= 1
                raise
        try:
            return self._idle.get(timeout=BORROW_TIMEOUT)
        except queue.Empty:
            raise RuntimeError(
                "no database session became available within %gs;"
                " every one of the %d in the pool is still busy"
                % (BORROW_TIMEOUT, self._maxsize))

    def release(self, session, discard=False):
        """
        Give a session back. discard=True throws it away instead, for
        when it may have been left in an unusable state.
        """
        if session is None:
            return
        if discard or self._closed:
            session.close()
            with self._lock:
                self._created -= 1
            return
        self._idle.put(session)

    def close(self):
        self._closed = True
        while True:
            try:
                self._idle.get_nowait().close()
            except queue.Empty:
                return


if __name__ == "__main__":
    import doctest
    doctest.testmod()
