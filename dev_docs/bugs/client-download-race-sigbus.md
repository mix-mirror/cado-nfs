# Clients sharing a working directory on NFS kill each other with SIGBUS

## Summary

When several `cado-nfs-client.py` processes share one working directory
on a network filesystem, the first program of each phase
(`polyselect`, then `las`) is downloaded by all of them at once into
the same path. Each downloads to a unique temporary name and then
renames it over the shared target. A client that is already *executing*
that target when another client's rename lands dies with **SIGBUS**,
and the workunit is reported as failed with error code `-7`.

The failures are concentrated at the start of each phase -- exactly
when the program for it is downloaded for the first time -- and then
stop, because later clients find the file already present and skip the
download. That makes them look like a transient glitch rather than a
race.

## Observed

A c130 factorization on 8 Grid'5000 nodes (59 clients, one shared
`tasks.workdir` on NFS) produced 10 failed workunits, all with error
code `-7`, all in the first moments of a phase, spread over seven
different machines:

    c130_polyselect1_1260-1890   -7   grvingt-31+4
    c130_polyselect1_1890-2520   -7   grvingt-16+7
    c130_polyselect1_2520-3150   -7   grvingt-31+7
    c130_polyselect1_3150-3780   -7   grvingt-9+3
    c130_polyselect1_3780-4410   -7   grvingt-12+8
    c130_sieving_750000-760000   -7   grvingt-14+8
    c130_sieving_760000-770000   -7   grvingt-15+4
    c130_sieving_790000-800000   -7   grvingt-12
    c130_sieving_800000-810000   -7   grvingt-16+3
    c130_sieving_810000-820000   -7   grvingt-15+3

The client log shows the download immediately before the death, with
an empty stderr:

    Downloading https://.../file/polyselect to .../client/download/polyselect165910178
    Setting executable flag for .../client/download/polyselect
    Running .../client/download/polyselect -P 65000 -N ... -admin 1260 -admax 1890 ...
    [...] Subprocess has PID 23186
    Command (...) resulted in exit code -7
    Child stderr:

and the next workunit on the same client runs the same program without
trouble:

    .../client/download/polyselect already exists, not downloading

## Reproducer

Deterministic, two nodes sharing an NFS home, no cado-nfs machinery
beyond one binary big enough to still be demand-paging:

```sh
D=$HOME/nfs-shared/bustest
BIN=<builddir>/polyselect/polyselect
N=7002225368837537027832023596485831272688308982276007693449270454855683576491474363461784762012463080438406772046201028336245352403

rm -rf "$D"; mkdir -p "$D"
cp "$BIN" "$D/polyselect"
cp "$BIN" "$D/polyselect.new"
sync

# node B: start a long run from the shared copy
ssh nodeB "$D/polyselect -P 65000 -N $N -degree 5 -t 2 \
           -admin 1000000 -admax 1200000 -incr 210 -nq 78125; \
           echo \$? > $D/rc" &
sleep 6

# node A: replace the executable underneath it
mv -f "$D/polyselect.new" "$D/polyselect"
wait
cat "$D/rc"
```

Result:

    Bus error
    135

135 = 128 + 7, i.e. SIGBUS -- the same signal the clients reported as
`-7`. On a local filesystem the same sequence is harmless: the running
process keeps the old inode. Across two NFS clients it is not, because
the pages of the executable are demand-paged from the server and the
file they came from has been replaced.

## Where it comes from

`cado-nfs-client.py`, in the download path:

```python
        if dlpath_tmp is not None:
            # We can't atomically rename-unless-dst-does-not-exist-yet.
            os.rename(dlpath_tmp, dlpath)
```

The comment is aware that the rename is unconditional. The existence
check that would avoid it happens much earlier ("%s already exists, not
downloading"), so between that check and this rename every client that
started at the same time decides to download, and each one in turn
replaces the file the others may be running.

## Fixed

Option 1 below, in this branch. A c100 on the same eight nodes, same
shared working directory, 64 clients, 293 workunits: no failures at
all, where the c130 above had five in polyselect alone.

## Suggested fixes

1. `os.link(dlpath_tmp, dlpath)` instead of the rename, then unlink the
   temporary. `link()` fails with `EEXIST` if the target exists, which
   is precisely the "atomically rename-unless-dst-does-not-exist-yet"
   the comment says is unavailable -- it is available, just not spelled
   `rename`. A client that loses the race keeps the copy that is
   already there, which is the same file.
2. Failing that, download to a path that is private to the client, so
   that no two clients ever write the same target. This costs one copy
   of each program per client.
3. Re-checking existence just before the rename narrows the window but
   does not close it.

Option 1 also fixes the case of a *stale* program being replaced
mid-run, which is the same race with a different trigger.

## Impact

Every failed workunit is resubmitted, so a run recovers by itself and
the only visible trace is a handful of `-7` failures in the first
minute of a phase. It does count against `tasks.maxfailed`, and on a
pool large enough that many clients start a phase simultaneously the
count grows with the number of clients rather than staying constant.
