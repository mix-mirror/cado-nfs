# Progress can exceed 100% because failed workunits are counted twice

## Summary

`Polysel1Task.get_achievement()` divides `wu_range_received` by the
`ad` span. `wu_range_received` counts the range of every workunit that
came back, **including the ones that came back as errors**. A failed
workunit is then resubmitted, and when the retry succeeds its range is
counted a second time. Progress therefore drifts above 1 in proportion
to the number of failures.

The code already knows:

```python
    def get_achievement(self):
        # Note that wu_range_received (like wu_received, by the way)
        # counts ERROR'd workunits as well !
        adspan = self.params["admax"] - self.params["admin"]
        return self.state["wu_range_received"] / adspan
```

## Observed

A c130 on 8 nodes, with `tasks.polyselect.admin=1260`,
`tasks.polyselect.admax=38e3`, so a span of 36740:

    wu_range_received    39890
    wu_failed            5
    achievement          1.0857376156777354

39890 / 36740 = 1.08574, which is the reported figure. The excess is
39890 - 36740 = 3150 = 5 x 630, i.e. exactly the five failed workunits
of range 630 each, counted once when they failed and once when their
retries succeeded.

(The five failures were themselves
[the NFS download race](client-download-race-sigbus.md), but any
source of failures produces this.)

## Reproducer

Any run with a failing polyselect workunit will do; the arithmetic is
deterministic given the failure. To force one without waiting for a
real fault, make a workunit fail by hand -- for instance run with a
`tasks.execpath` whose `polyselect` exits non-zero on the first call
-- and compare `wu_range_received` against
`admax - admin` in the database afterwards:

    sqlite3 <name>.db 'SELECT kkey,value FROM polyselect1
                       WHERE kkey IN ("wu_range_received","wu_failed")'

`wu_range_received` will exceed the span by the range of each failed
workunit.

## Suggested fix

Do not add the range of a workunit that came back as an error, or
subtract it again when the retry is counted. The same `while here`
note applies to `wu_received`, which the comment mentions and which
feeds the "N / M workunits back" figure.

`SievingTask` and `PolyselJLTask` compute achievement the same way and
are presumably affected in the same proportion.

## Note

The dashboard currently clamps the *displayed* percentage to 100% and
keeps the true figure in the tooltip, since a ring that is visibly
full must not be labelled 108.6%. That is a presentation patch over
this bug, not a fix for it.
