#!/usr/bin/env bash

# Pin the output of the deterministic half of the strategies chain.
#
# gst, in all four of its modes, and finalst draw no random numbers: they
# are pure functions of their input files. This test runs them on fixed
# inputs and compares the result against references committed beside it,
# so that reworking sieve/strategies has a net that catches any change of
# behaviour rather than only a change of exit status.
#
# The other half of the chain -- gfm and benchfm -- measures wall-clock
# time and cannot be pinned this way; test_strategies_pinned_notime.sh
# covers it modulo the timing columns.
#
# After a deliberate change of behaviour, re-run with REGEN=1 to rewrite
# the references in the source tree, and read the diff before committing.

set -e

: ${CADO_NFS_SOURCE_DIR:?missing}
: ${CADO_NFS_BINARY_DIR:?missing}
: ${wdir:?missing}

GST="$CADO_NFS_BINARY_DIR/sieve/strategies/gst"
FINALST="$CADO_NFS_BINARY_DIR/sieve/strategies/finalst"
here="$CADO_NFS_SOURCE_DIR/tests/sieve/strategies"
ref="$here/expected"
sample="$here/sample_methods.txt"

lim=1024
lpb=16
mfb=26
ncurves=3

# fbb = ceil(log2(lim+1)) = 11 here, so gst treats a cofactor of fewer
# than 2*fbb-1 = 21 bits as prime and hands it a trivial strategy. The
# band [21, mfb] is where the generator actually does work, and it is
# where the references below are taken.
rmin=21

# The pairs pinned verbatim: the two ends of the diagonal, one cell off
# it, and one strongly asymmetric cell.
pinned_pairs="21_21 21_26 24_25 26_26"

cd "$wdir"
rc=0

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

label() { # concatenate files, each preceded by its basename
    local f
    for f in "$@" ; do
        echo "== $(basename "$f")"
        cat "$f"
    done
}

############################################################
echo "### gst -gdc"
mkdir -p dec
for r in $(seq $rmin $mfb) ; do
    $GST -gdc -lim0 $lim -mfb0 $r -out dec/decomp_${lim}_$r
done
label dec/decomp_${lim}_* > decomps.txt
compare decomps.txt decomps.txt

############################################################
echo "### gst -gst_r"
mkdir -p pre
for r in $(seq $rmin $mfb) ; do
    $GST -gst_r -lim0 $lim -lpb0 $lpb -r0 $r -ncurves $ncurves \
        -in "$sample" -decomp dec/decomp_${lim}_$r -out pre
done
# the trivial sizes below rmin are needed by the -gst step further down
for r in $(seq 0 $((rmin-1))) ; do
    $GST -gst_r -lim0 $lim -lpb0 $lpb -r0 $r -ncurves $ncurves \
        -in "$sample" -decomp /dev/null -out pre
done
label $(for r in $(seq $rmin $mfb) ; do echo pre/strategies${lim}_$r ; done) > oneside.txt
compare oneside.txt oneside.txt

############################################################
echo "### gst -gst (one pair at a time)"
mkdir -p pair
for r0 in $(seq $rmin $mfb) ; do
    for r1 in $(seq $rmin $mfb) ; do
        $GST -gst -lim0 $lim -lim1 $lim -r0 $r0 -r1 $r1 -in pre -out pair
    done
done
label $(for p in $pinned_pairs ; do echo pair/strategies_$p ; done) > pairs.txt
compare pairs.txt pairs.txt

############################################################
echo "### gst, default mode (whole matrix in one go)"
mkdir -p mat
$GST -lim0 $lim -lpb0 $lpb -mfb0 $mfb -lim1 $lim -lpb1 $lpb -mfb1 $mfb \
    -ncurves $ncurves -in "$sample" -decomp dec -out mat > /dev/null

n=$(ls mat | wc -l)
if [ "$n" != $(( (mfb+1) * (mfb+1) )) ] ; then
    echo "gst default mode wrote $n files, expected $(( (mfb+1)*(mfb+1) ))" >&2
    rc=1
fi
label $(for p in $pinned_pairs ; do echo mat/strategies_$p ; done) > matrix.txt
compare matrix.txt matrix.txt

# the whole matrix, as one digest: too large to pin verbatim, but a
# change anywhere in it should still be noticed.
( cd mat && cat $(ls | sort) ) | sha256sum | cut -d' ' -f1 > matrix.sha256
compare matrix.sha256 matrix.sha256

############################################################
echo "### the two routes to a pair must agree on the method chains"
# They are not byte-identical: -gst reads the one-sided strategies back
# from a file that stores them with six decimals, so the probabilities
# and times it derives differ from the in-memory ones in the last digit.
# The chains of methods themselves must match exactly.
for r0 in $(seq $rmin $mfb) ; do
    for r1 in $(seq $rmin $mfb) ; do
        a=pair/strategies_${r0}_${r1}
        b=mat/strategies_${r0}_${r1}
        if ! diff -q <(grep -v '^Probability:\|^Time:' $a) \
                     <(grep -v '^Probability:\|^Time:' $b) > /dev/null ; then
            echo "cell ($r0,$r1): -gst and the default mode disagree on the methods" >&2
            rc=1
        fi
    done
done

############################################################
echo "### finalst"
# A real cofactor distribution is sparse -- most size pairs never occur --
# so the reference stays small without losing any of the selection logic.
for r0 in 21 23 24 26 ; do
    for r1 in 21 23 24 26 ; do
        echo "$r0 $r1 $(( 1 + (r0 * 7 + r1 * 3) % 11 )) 1"
    done
done > cofactors.stats

# two values of -t, to land on two different points of the slope
# selection in compute_best_strategy().
for t in 0.01 100 ; do
    $FINALST -st mat -dist cofactors.stats -t $t -mfb0 $mfb -mfb1 $mfb \
        -out final_st_$t > /dev/null
    compare final_st_$t final_st_$t
done

exit $rc
