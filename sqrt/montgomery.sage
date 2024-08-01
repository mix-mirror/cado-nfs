
%load_ext autoreload
%autoreload 2

import subprocess
import cado_sage
from cado_sage import CadoPolyFile
from cado_sage import CadoPurgedFile
from cado_sage import CadoIndexFile
from cado_sage import CadoIdealsMapFile
from cado_sage import CadoIdealsDebugFile
from cado_sage import CadoNumberTheory
from cado_sage.bwc_sage import *

cado_sage.set_verbose(True)

exponent=65537
wdir="/tmp/blah"
name="c30"
number_of_solutions=3
par = BwcParameters(m=8,n=4,p=exponent,wordsize=64)


matrixfile = f"{wdir}/{name}.sparse.bin"
solution_files = [f"{wdir}/{name}.bwc/K.sols{i}-{i+1}.0.txt"
                  for i in range(number_of_solutions)]
purgedfile = f"{wdir}/{name}.purged.gz"
indexfile = f"{wdir}/{name}.index.gz"
polyfile = f"{wdir}/{name}.poly"
idealsmapfile = f"{wdir}/{name}.ideals.gz"
# This is created by explain_indexed_relation -python -all
idealsdebugfile = f"{wdir}/{name}.debug-ideals.txt"


# take a centered representative
centered = lambda x, p: (x-p) if 2*x > p else x

M = BwcMatrix(par, matrix=matrixfile, wdir=wdir)
M.read()
A =  M.M.dense_matrix()
ker = matrix([[ZZ(c) for c in open(f).readlines()] for f in solution_files])
assert matrix(GF(par.p), ker) * A == 0


abpairs = CadoPurgedFile(purgedfile)
abpairs.read_abpairs()


relsets = CadoIndexFile(indexfile)
relsets.read()
R = relsets.matrix()

poly = CadoPolyFile(polyfile)
poly.read()

K0 = poly.K[0]; alpha0 = K0.gen()
K1 = poly.K[1]; alpha1 = K1.gen()

unit_rank = [sum(K.signature()) - 1 for K in poly.K]

if sum(unit_rank) > number_of_solutions:
    raise ValueError(f"number_of_solutions={number_of_solutions} will not" +
                     f" disambiguate units (ranks: {unit_rank})")


nt = CadoNumberTheory(poly)

# It's perhaps a bit too much to ask, because sagemath will factor the
# discriminant, here.
OK0, OK1 = nt.maximal_orders()
J0, J1 = nt.J()




L1 = nt[1].LogDriftMap()

E = matrix([ L1(a-b*alpha1) for a,b,_ in abpairs ])


# all the valuations that we have for our current e-th power.

ideal_valuations_matrix = matrix(ZZ,ker*M.M.change_ring(ZZ) / exponent)
log_embeddings_matrix = matrix(ZZ, ker * R * E / exponent * 2^64)

vals = block_matrix([[ideal_valuations_matrix, log_embeddings_matrix]])


ideals_map = CadoIdealsMapFile(idealsmapfile)
ideals_map.read()


# read the ideals in sage format as well. This takes about 4 seconds for
# a {name}...
all_ideals = CadoIdealsDebugFile(poly, idealsdebugfile)
all_ideals.read()


schirokauer_block = matrix([vector(sum([list(f(a-b*alpha1)) for f in m],[])) for a,b,rel in abpairs])

assert (ker * R * schirokauer_block).rank() <= sum(unit_rank)

# This is such that the product of (a,b) pairs with these exponents is
# now guaranteed to be an ell-th power.
ker_power = S.left_kernel().basis()[0] * ker
ker_power = vector([centered(ZZ(x), par.p) for x in list(ker_power)])

# and this is the drift we're going to compete against.
drift = ker_power * R * E
