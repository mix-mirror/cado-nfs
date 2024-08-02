# You probably want to run this with export PYTHONPATH=./tests/sagemath

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
from cado_sage import CadoMontgomeryReductionProcess
from cado_sage.bwc import BwcParameters, BwcMatrix

cado_sage.set_verbose(True)

exponent=23
wdir="/tmp/blah"
name="c30"
number_of_solutions=3
par = BwcParameters(m=8,n=4,p=exponent,wordsize=64)

# not sure we'll ever need more.
prec = 53


matrixfile = f"{wdir}/{name}.sparse.bin"
solution_files = [f"{wdir}/{name}.bwc.{exponent}/K.sols{i}-{i+1}.0.txt"
                  for i in range(number_of_solutions)]
purgedfile = f"{wdir}/{name}.purged.gz"
indexfile = f"{wdir}/{name}.index.gz"
polyfile = f"{wdir}/{name}.poly"
idealsmapfile = f"{wdir}/{name}.ideals.gz"
# This is created by explain_indexed_relation -python -all
idealsdebugfile = f"{wdir}/{name}.debug-ideals.txt"


M = BwcMatrix(par, matrix=matrixfile, wdir=wdir)
M.read()
ker = matrix([[ZZ(c) for c in open(f).readlines()] for f in solution_files])
assert ker * M.M == 0


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


ideals_map = CadoIdealsMapFile(idealsmapfile)
ideals_map.read()


# read the ideals in sage format as well. This takes about 4 seconds for
# a {name}...
all_ideals = CadoIdealsDebugFile(poly, idealsdebugfile)
all_ideals.read()


# We're only really interested in the valuations per number field. Since
# we now have all the information available, we can split the matrix in
# two blocks, one per number field.

matrices_per_field = []
ideals_per_field = []
for K in poly.K:
    indices = [i for i,j in enumerate(ideals_map)
               if all_ideals[j].number_field() == K]
    MK = matrix(M.M.ncols(), len(indices),
                entries={(j,i):1 for i,j in enumerate(indices)})
    matrices_per_field.append((M.M*MK).lift_centered())
    ideals_per_field.append([ all_ideals[ideals_map[i]] for i in indices ])
split_valuations_matrix = block_matrix(1,2,matrices_per_field)



# Use Schirokauer maps to guarantee an ell-th power. We do this modulo
# both fields, but it can well be that it's useless modulo one of them.
# It's not much of an issue, since in that case we're going to have a
# zero-rank block and that's all

# we need to precompute the maps and not create them multiple times.
sm = [ n.schirokauer_maps(exponent) for n in nt ]
def sm_apply(m, x):
    L = [f(x) for f in m]
    return flatten(L, ltypes=(type(L[0]),))

schirokauer_block = block_matrix([[matrix([sm_apply(sm[i],a - b * n.K.gen())
                                           for a,b,rel in abpairs])
                                   for i,n in enumerate(nt)]])
S = ker * R * schirokauer_block
assert S.subdivision(0,0).rank() <= unit_rank[0]
assert S.subdivision(0,1).rank() <= unit_rank[1]

if False:
    # The code here is much simpler, but it is also a lot slower. The
    # good thing is that we do compute the same thing!
    sm2 = [n.schirokauer_map_simple(exponent) for n in nt]
    schirokauer_block2 = matrix([vector(sum([list(sm2[i](a - b * n.K.gen()))
                                        for i,n in enumerate(nt)
                                        ],[]))
                                 for a,b,rel in abpairs])
    S2 = ker * R * schirokauer_block2

    assert S.left_kernel() == S2.left_kernel()


# This is such that the product of (a,b) pairs with these exponents is
# now guaranteed to be an ell-th power. It's important to take a centered
# lift, as it minimizes the valuations that we get later on!
ker_power = (S.left_kernel().basis()[0] * ker).lift_centered()

# Now let's focus on one side in particular.
side = 1

K = poly.K[side]
OO = nt[side].OrderCastMap()
OK = K.maximal_order()

ideals = ideals_per_field[side]

# all the valuations that we have for our current e-th power.
valuations = ker_power * split_valuations_matrix.subdivision(0,side) / exponent

L = nt[side].LogMap(prec)

E = matrix([ L(a - b * K.gen()) for a,b,rel in abpairs ])

# these are the log embeddings of our root
embeddings = ker_power * R * E / exponent

MM = CadoMontgomeryReductionProcess(poly, side, ideals, valuations, embeddings)
MM.status()

while (MM.nplus + MM.nminus)[1] > 100:
    MM.one_reduction_step(256)
    MM.status()

gamma = MM.accumulated.prod()

# even for exponent=23, this takes an absurdly long time. In production,
# we'll never do this.
print("computing the absurdly big algebraic number (LONG!)")
big_power = prod([ (ab[0] - ab[1] * alpha1)**(ker_power * R)[i] for i,ab in enumerate(abpairs) ])

print("computing the remaining bit after all cancellations")
rest = big_power / gamma**exponent

print("computing the final e-th root")
delta = rest.nth_root(exponent)

# then (gamma * delta) is an e-th root of big_power
print("Final check: ", big_power / (gamma*delta)**exponent)
