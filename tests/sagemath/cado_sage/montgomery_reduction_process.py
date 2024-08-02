from sage.modules.free_module_element import vector
from sage.functions.generalized import sign
from sage.functions.other import ceil
from sage.misc.functional import round
from sage.structure.factorization import Factorization
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix
from sage.functions.log import exp
from sage.matrix.special import vandermonde
from sage.rings.complex_double import CDF
from sage.rings.complex_mpfr import ComplexField
import math


class CadoMontgomeryReductionProcess(object):
    def __init__(self, poly, side,
                 ideals, valuations, embeddings):
        """
        with respect to the given side, and the corresponding matrices
        that give the valuations and log embeddings for that field,
        return an algebraic number that matches these as best as we can.
        """

        self.current = { I : valuations[i] for i, I in enumerate(ideals)
                        if valuations[i] }
        self.accumulated = 1
        self.embeddings = embeddings
        self.nplus = vector([0, 0])
        self.nminus = vector([0, 0])
        self.K = poly.K[side]
        self.nt = poly.nt[side]
        self.OK = self.K.maximal_order()
        self.L = poly.nt[side].LogMap()
        self.poly = poly

    def skewed_LLL(self, M, skew, algorithm='pari'):
        """
        returns an LLL basis where the coefficients of each vector are
        skewed according to the skewness of the number field polynomial

        TODO: what if the skew matrix is a complex matrix? I guess that
        we would have to use the conjugate_transpose, right? I'm slightly
        puzzled that this implies losing the nice structure of a Hankel
        matrix with the Newton sums when we skew by a Vandermonde matrix.

        TODO: we should be accepting an inner product matrix of the
        ambient space, not just a vector of weights. In order to make it
        work with fplll, this requires the computation of the Cholesky
        decomposition.
        """

        D = skew

        # it seems that sage won't let us do an LLL reduction just based
        # on a Gram matrix while using the underlying fplll code.
        # so if I insist on using fpLLL, what I need to do is
        if algorithm == 'fplll':
            # If the scaling weights are below one, we must pay attention
            # to what we're doing!
            dm = min(D.diagonal())
            if dm < 1:
                D = D / dm
            D = D.apply_map(round)
            return matrix(M.base_ring(), (M * D).LLL() * D ** -1)

        # sage does expose an LLL_gram, but it's from pari, and it
        # behaves quite oddly.
        elif algorithm == 'pari':
            D2 = D * D.transpose()
            for i in range(20):
                # It's really abominable, but apparently we _have_ to do
                # that and recompute the gram matrix several times,
                # otherwise pari outputs bogus stuff. Oh dear.
                G = (M * D2 * M.transpose())
                T = G.LLL_gram()
                M = T.transpose() * M
                if T.norm() < 2:
                    break
            return M

    def pick_an_ideal_product_to_kill(self, num_or_den=1, norm_max_bits=64):
        """
        return an ideal (in factored form) that we will attempt to kill,
        with a bound on the number of bits of its norm
        """
        cumulative_norm_bits = 0
        # i = 0
        # V = vector([0]*len(ideals))
        pool = []
        norm_desc = []
        for I, emax in reversed(self.current.items()):
            if cumulative_norm_bits >= norm_max_bits:
                break
            if sign(emax) != num_or_den:
                continue
            nn = math.log(I.norm(), 2)
            cap = ceil((norm_max_bits - cumulative_norm_bits) / nn)
            e = min(emax * num_or_den, cap)
            pool.append((I, e))
            # V[i] += e
            cumulative_norm_bits += e * nn
            norm_desc.append(f"{nn:.0f}" + (f"*{e}" if e > 1 else ""))
        norm_desc = " + ".join(norm_desc)
        print(f"Picked ideal (sign {num_or_den})"
              + f" has {cumulative_norm_bits:.2f}-bit norm ({norm_desc})")
        return Factorization(pool)

    def status(self):
        """
        Some stats. Note that this incidentally updates nplus and nminus
        """
        offset_bits = self.embeddings.norm() / math.log(2)
        self.nplus = vector([0, 0])
        self.nminus = vector([0, 0])
        for I, v in self.current.items():
            if v > 0:
                self.nplus += vector([1, v * math.log(I.norm(), 2)])
            else:
                self.nminus += vector([1, -v * math.log(I.norm(), 2)])

        print(f"offset bits {offset_bits:.2f},"
              + f" num bits {self.nplus[1]:.2f} ({self.nplus[0]:.0f} ideals),"
              + f" denom bits {self.nminus[1]:.2f} ({self.nminus[0]:.0f} ideals)")

    def accumulate(self, g, num_or_den, hint=Factorization([])):
        """
        Take action, and register a new term for the product, updating
        the current state accordingly
        """
        self.accumulated *= Factorization([(g, num_or_den)])
        print(f"Picked gen has {math.log(g.norm().abs(), 2):.2f}-bit norm")
        self.embeddings -= num_or_den * self.L(g)

        # TODO: we should not be relying on the maximal order here.
        discovered = (self.OK.fractional_ideal(g) / hint.prod()).factor()
        ideal_factorization = hint * discovered

        for I, e in ideal_factorization:
            # Ibits = math.log(I.norm(), 2)
            if I not in self.current:
                # print(f"Just added a {Ibits:.2f}-bit ideal")
                self.current[I] = 0
            self.current[I] -= num_or_den * e
            if self.current[I] == 0:
                # print(f"Just got rid of {Ibits:.2f}-bit ideal")
                del self.current[I]

    def one_reduction_step(self, bits=64):
        num_or_den = 1 if self.nplus[1] > self.nminus[1] else -1

        hint = self.pick_an_ideal_product_to_kill(num_or_den, bits)
        I = hint.prod()
        assert I.is_integral()

        # by doing LLL here, we tolerate the norm of each generator to
        # grow by as much as a constant C_K, which we can compute.
        gens0 = I.basis()
        M0 = matrix([list(c) for c in gens0])

        # There does seem to be a computational advantage in doing the
        # reduction in two steps.
        d = self.K.degree()
        s = self.poly.skewness
        skew0 = diagonal_matrix([(s ** (i - (d - 1)/2)) for i in range(d)])
        M1 = self.skewed_LLL(M0, skew0)

        gens1 = [self.K(list(g)) for g in M1]

        # The theory is that the skewed norm of the vector M1[0] is
        # bounded as follows
        assert float((M1[0]*skew0).norm()) < float(2**((d-1)/4)*M0.determinant()**(1/d))


        # We now have alternative, somewhat smaller generators of our
        # ideal. They're finely skewed, which contributes to making the
        # norm a bit smaller than if we had neglected that aspect.

        # Next we want the embeddings to be small, which we translate
        # into the requirement that the vector that we end up with is as
        # close as we can to the same skewed hypersphere as the one the
        # current embeddings of the target lie on.

        modules = self.nt.modules_of_embeddings_from_log_embeddings(self.embeddings)

        f = self.K.defining_polynomial()
        prec = 53

        if prec == 53:
            C = CDF
        else:
            C = ComplexField(prec)


        V = vandermonde(f.roots(C, multiplicities=False)).transpose()

        # note that V * V.transpose() is the Hankel matrix whose
        # coefficients are the Newton sums, which we can compute fairly
        # easily via the expansion of (xf'/f) in powers of 1/x

        # alas, we're rather interested by the comparison with the
        # hypersphere, and this kills the nice properties of these
        # expression. The polynomial that gives the target embeddings is
        # unknown anyways. We have to consider a skewed V:

        # if we take:
        # W = V * diagonal_matrix(self.nt.real_representation_of_embeddings(g))**-1
        # then W is so that vector(g.list()) * W is the all-ones vector

        # but of course, the W that we want is the one that compares to
        # the target number!
        W = V * diagonal_matrix(modules)**-num_or_den

        # So. At this point, what we have to do is to give a skewed LLL
        # reduction with respect to this matrix W
        M2 = self.skewed_LLL(M1, skew=W)
        gens2 = [self.K(list(g)) for g in M2]

        # Just pick the smallest one.
        # Now we're still very visibly doing something stupid, which is
        # reducing the same lattice consecutively with two different
        # norms. We should rather take both into consideration.
        g = gens2[0]

        drift_bits = (vector(g.list()) * W).norm().log(2)
        print(f"Picked gen has {drift_bits:.2f}-bit drift")

        self.accumulate(g, num_or_den, hint=hint)

    def reduction_process(self):
        pass
