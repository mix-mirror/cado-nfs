#!/usr/bin/env bash

# Pin what can be pinned of gfm and benchfm.
#
# These two measure wall-clock time and write it into their output, so
# they never produce the same bytes twice. With -seed fixed, though,
# everything *else* they write is reproducible: the methods they select
# and the probabilities they measure. This test therefore compares their
# output with the timing column stripped off.
#
# Not every mode survives that treatment:
#
#   gfm (no -ch), benchfm -p -t     reproducible once the timings are
#                                   stripped -- pinned here.
#   benchfm -p, benchfm -f, gfm     no measurement involved, or only a
#   -fch                            pinned file as input: reproducible as
#                                   they stand -- pinned verbatim here.
#   gfm -ch                         NOT reproducible, even stripped: the
#                                   convex hull ranks methods by
#                                   (probability, time), so timing noise
#                                   changes which methods survive. Eight
#                                   runs at one seed gave three different
#                                   answers, with hulls of different
#                                   sizes. Only smoke-checked below; the
#                                   same hull code is pinned properly
#                                   through -fch, which reads a fixed file.
#
# Re-run with REGEN=1 to rewrite the references after a deliberate change.

set -e

: ${CADO_NFS_SOURCE_DIR:?missing}
: ${CADO_NFS_BINARY_DIR:?missing}
: ${wdir:?missing}

GFM="$CADO_NFS_BINARY_DIR/sieve/strategies/gfm"
BENCHFM="$CADO_NFS_BINARY_DIR/sieve/strategies/benchfm"
here="$CADO_NFS_SOURCE_DIR/tests/sieve/strategies"
ref="$here/expected"
sample="$here/sample_methods.txt"

# The sieve region is chosen so that c stays well clear of B1: gfm aborts
# in stage2_make_plan() when B2 = c*B1 lands within 30 of B1*B1, which the
# default region does reach.
region="-b1min 100 -b1max 400 -b1step 100 -cmin 10 -cmax 40 -cstep 10"

cd "$wdir"
rc=0

# drop the third |-delimited field, which holds the measured timings
strip_timings() { awk -F'|' '{print $1 "|" $2 "|"}' "$1" ; }

compare() { # <name> <actual>
    local name="$1" actual="$2"
    if [ "$REGEN" ] ; then
        mkdir -p "$ref"
        cp "$actual" "$ref/$name"
        echo "REGEN: wrote $ref/$name"
        return
    fi
    if [ ! -f "$ref/$name" ] ; then
        echo "missing reference $ref/$name (re-run with REGEN=1)" >&2
        rc=1
        return
    fi
    if ! diff -u "$ref/$name" "$actual" > "diff_$name.txt" ; then
        echo "$name: output differs from the pinned reference" >&2
        head -40 "diff_$name.txt" >&2
        rc=1
    fi
}

############################################################
echo "### gfm, one family at a time, timings stripped"
for m in PM1 PP1-27 ECM-M12 ECM-M16 ECM-B12 ; do
    $GFM -lb 16 -ub 19 -m $m $region -seed 1 -out gfm_$m
    strip_timings gfm_$m > gfm_$m.notime
    compare gfm_$m.notime gfm_$m.notime
done

############################################################
echo "### gfm -fch (convex hull of a pinned file)"
$GFM -fch -fch_in "$sample" -fch_out gfm_fch.txt > /dev/null
compare gfm_fch.txt gfm_fch.txt

############################################################
echo "### gfm -ch (smoke only, see the header)"
$GFM -lb 16 -ub 18 -m ECM-M12 $region -seed 1 -ch -out gfm_ch.txt
if [ ! -s gfm_ch.txt ] ; then
    echo "gfm -ch produced an empty hull" >&2
    rc=1
fi
# whatever survives must be one of the methods that were swept
$GFM -lb 16 -ub 18 -m ECM-M12 $region -seed 1 -out gfm_noch.txt
while read -r m c b1 b2 rest ; do
    if ! cut -d'|' -f1 gfm_noch.txt | grep -qx "$m $c $b1 $b2 " ; then
        echo "gfm -ch returned a method that was never swept: $m $c $b1 $b2" >&2
        rc=1
    fi
done < <(cut -d'|' -f1 gfm_ch.txt)

############################################################
echo "### benchfm -p (probabilities only: no measurement in the output)"
$BENCHFM -in "$sample" -p -lb 10 -N 200 -seed 2 -out benchfm_p.txt
compare benchfm_p.txt benchfm_p.txt

############################################################
echo "### benchfm -p -t, timings stripped"
$BENCHFM -in "$sample" -p -t -lb 10 -N 200 -seed 2 -out benchfm_pt
strip_timings benchfm_pt > benchfm_pt.notime
compare benchfm_pt.notime benchfm_pt.notime

# -p alone and -p -t must agree on the probabilities they measure: same
# seed, same input, and the timing pass does not touch them.
if ! diff -q <(strip_timings benchfm_p.txt) benchfm_pt.notime > /dev/null ; then
    echo "benchfm -p and benchfm -p -t disagree on the probabilities" >&2
    diff <(strip_timings benchfm_p.txt) benchfm_pt.notime | head -10 >&2
    rc=1
fi

############################################################
echo "### benchfm -f (filtering a pinned file)"
$BENCHFM -in "$sample" -f 3 -lb 10 -out benchfm_f3.txt
compare benchfm_f3.txt benchfm_f3.txt

exit $rc
