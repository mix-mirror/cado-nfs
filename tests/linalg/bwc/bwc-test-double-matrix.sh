#!/usr/bin/env bash

set -e

tmp=$(mktemp -d /tmp/cado-nfs.XXXXXXXXXXXXXX)

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

TESTCASE="123 251 172 STRETCH STRETCH"  # for thr=2x2
TESTCASE="159 193 129 STRETCH STRETCH"  # for thr=1x1
bindir="$CMAKE_CURRENT_BINARY_DIR"
thr=2x2
seed=1

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
    else
        usage
    fi
done

# generate the matrices
`dirname "$0"`/create-fake-dble-matrices.pl --save-matrices-list $tmp/matrices.txt --dstdir $tmp --seed $seed $TESTCASE

# form the comma-separated list of matrices
matrices=""
set -- `cat $tmp/matrices.txt`
for m in "$@" ; do
    if [ "$matrices" ] ; then matrices="$matrices," ; fi
    matrices="$matrices$m"
done

common=(
        multi_matrix=1
        matrix=$matrices
        wdir=$tmp
        mn=64
        balancing_options=reorder=none,skip_decorrelating_permutation=1
        thr=$thr
    )
$bindir/prep "${common[@]}" seed=$seed
$bindir/secure "${common[@]}" interval=1

