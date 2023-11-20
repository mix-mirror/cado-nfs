#!/usr/bin/env bash

set -e

tmp=$(mktemp -d /tmp/cado-nfs.bwc-test-double-matrix.XXXXXXXXXXXXXX)

if ! [ "$CADO_DEBUG" ] ; then
    trap "rm -rf $tmp" EXIT
else
    set -x
fi

usage() {
    echo "Options:" >&2
    echo "  --bindir <path>     places where bwc binaries are found" >&2
    echo "  --test <argument string>    arguments passed to script that generates the test case" >&2
    exit 1
}

# TESTCASE="123 251 172 STRETCH STRETCH"  # for thr=2x2
TESTCASE="159 193 129 STRETCH STRETCH"  # for thr=1x1
bindir="$CMAKE_CURRENT_BINARY_DIR"
thr=2x2
seed=1

prime=2
m=64
n=64

while [ $# -gt 0 ] ; do
    if [ "$1" = "--bindir" ] ; then
        shift
        bindir="$1"
        shift
    elif [ "$1" = "--test" ] ; then
        shift
        TESTCASE="$1"
        shift
    elif [ "$1" = "--seed" ] ; then
        shift
        seed="$1"
        shift
    elif [ "$1" = "--thr" ] ; then
        shift
        thr="$1"
        shift
    elif [ "$1" = "--prime" ] ; then
        shift
        prime="$1"
        shift
    elif [ "$1" = "--m" ] ; then
        shift
        m="$1"
        shift
    elif [ "$1" = "--n" ] ; then
        shift
        n="$1"
        shift
    else
        usage
    fi
done

# generate the matrices
`dirname "$0"`/create-fake-dble-matrices.pl --prime $prime --save-matrices-list $tmp/matrices.txt --dstdir $tmp --seed $seed $TESTCASE

# form the comma-separated list of matrices
matrices=""
set -- `cat $tmp/matrices.txt`
for mat in "$@" ; do
    if [ "$matrices" ] ; then matrices="$matrices," ; fi
    matrices="$matrices$mat"
done

pmn=( prime=$prime m=$m n=$n)

case "$prime" in
    2) lingen=lingen_b64;;
    *) lingen=lingen_p1;;
esac

eval "${pmn[@]}"

common=(
        multi_matrix=1
        matrix=$matrices
        "${pmn[@]}"
        wdir=$tmp
        balancing_options=reorder=none,skip_decorrelating_permutation=1
        thr=$thr
    )
$bindir/prep "${common[@]}" seed=$seed
$bindir/secure "${common[@]}"
$bindir/krylov "${common[@]}" ys=0..$n
$bindir/acollect wdir=$tmp "${pmn[@]}" remove_old=1
afiles=(`find $tmp -name "A0-*.0-*"`)
if [ "${#afiles[@]}" != 1 ] ; then
    echo "Error: want exactly one A file for lingen" >&2
    exit 1
fi
afile="${afiles[0]}"
ffile="$afile.gen"
$bindir/$lingen "${pmn[@]}" wdir=$tmp afile="$afile" ffile=F rhs=none split-output-file=1
$bindir/mksol "${common[@]}" solutions=0-$n
$bindir/gather "${common[@]}" solutions=0-$n
