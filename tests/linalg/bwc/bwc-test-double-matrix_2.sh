#!/bin/bash
# https://gitlab.inria.fr/cado-nfs/cado-nfs/-/issues/30070
p=`pwd`
d=`mktemp -d`
cd $d
wget https://gitlab.inria.fr/cado-nfs/cado-nfs/uploads/d25129cf3974d86e7a0fd10fd71ced1b/R.rw.bin
wget https://gitlab.inria.fr/cado-nfs/cado-nfs/uploads/31d04463ed3d8b281bb199d35a96928d/R.dense.rw.bin
wget https://gitlab.inria.fr/cado-nfs/cado-nfs/uploads/ae1220c9cfff6e0fda2058bb29a610d1/R.dense.cw.bin
wget https://gitlab.inria.fr/cado-nfs/cado-nfs/uploads/fe42bdf5ed7111df16e4ff5c9636c797/R.dense.bin
wget https://gitlab.inria.fr/cado-nfs/cado-nfs/uploads/9cad98eb94b0b0f58c83929859ab2f9b/R.cw.bin
wget https://gitlab.inria.fr/cado-nfs/cado-nfs/uploads/844849ec8988228dd9ebc41b60ac86b5/R.bin
wget https://gitlab.inria.fr/cado-nfs/cado-nfs/uploads/d2325acbe43b04b0cbca7def5c5968d1/L.rw.bin
wget https://gitlab.inria.fr/cado-nfs/cado-nfs/uploads/76752bc10757f9f97483e42a6f5db363/L.cw.bin
wget https://gitlab.inria.fr/cado-nfs/cado-nfs/uploads/fd79194a888208fb667606878d5c429b/L.bin
echo $d
cd $p
./build/`hostname`/linalg/bwc/prep multi_matrix=1 matrix=$d/L.bin,$d/R.bin wdir=$d/bwc mn=64 balancing_options=reorder=columns thr=2x2
./build/`hostname`/linalg/bwc/secure multi_matrix=1 matrix=$d/L.bin,$d/R.bin wdir=$d/bwc mn=64 thr=1x1 prime=2 balancing_options=reorder=none,skip_decorrelating_permutation=1 
./build/`hostname`/linalg/bwc/krylov multi_matrix=1 matrix=$d/L.bin,$d/R.bin wdir=$d/bwc mn=64 thr=1x1 prime=2 balancing_options=reorder=none,skip_decorrelating_permutation=1 ys=0..64
./build/`hostname`/linalg/bwc/acollect wdir=$d/bwc mn=64 prime=2 ys=0..64 --remove-old
afiles=(`find $d/bwc -name "A0-*.0-*"`)
if [ "${#afiles[@]}" != 1 ] ; then
    echo "Error: want exactly one A file for lingen" >&2
    exit 1
fi
afile="${afiles[0]}"
ffile="$afile.gen"
./build/`hostname`/linalg/bwc/lingen_b64 wdir=$d/bwc mn=64 prime=2 afile="$afile" ffile=F rhs=none split-output-file=1
./build/`hostname`/linalg/bwc/mksol multi_matrix=1 matrix=$d/L.bin,$d/R.bin wdir=$d/bwc mn=64 thr=1x1 prime=2 balancing_options=reorder=none,skip_decorrelating_permutation=1 solutions=0-64
./build/`hostname`/linalg/bwc/gather multi_matrix=1 matrix=$d/L.bin,$d/R.bin wdir=$d/bwc mn=64 thr=1x1 prime=2 balancing_options=reorder=none,skip_decorrelating_permutation=1 solutions=0-64
