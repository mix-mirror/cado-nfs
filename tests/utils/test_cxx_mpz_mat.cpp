#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <utility>
#include <map>
#include <algorithm>

#include <gmp.h>

#include "gmp_aux.h"
#include "mpz_mat.h"

// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_gmp_randstate state;

    if (argc > 1) {
        unsigned long seed;
        seed = strtoul(argv[1], nullptr, 0);
        gmp_randseed_ui(state, seed);
    }

    std::map<unsigned long, cxx_mpz_mat> v;

    mpz_t det;
    mpz_init(det);
    unsigned long const p = 1009;

    /* This generates many matrices at random, and puts in the std::map
     * above the relationship with the determinant of the reduction mod
     * 1009 of the leading submatrix of their HNF
     */
    for(int i = 0 ; i < 10 ; i++) {
        cxx_mpz_mat M;

        mpz_mat_realloc(M, gmp_urandomb_ui(state, 4) + 2, gmp_urandomb_ui(state, 4) + 2);
        for(unsigned int i = 0 ; i < M->m ; i++) {
            for(unsigned int j = 0 ; j < M->n ; j++) {
                mpz_set_si(mpz_mat_entry(M, i, j), (static_cast<long>(gmp_urandomb_ui(state, 16)) - 32768));
            }
        }
        mpz_mat_mod_ui(M, M, p);
        cxx_mpz_mat M1 = M, T, M2;
        mpz_mat_hermite_form(M1, T);
        // mpz_mat_fprint(stdout, M1);
        unsigned int const d = std::min(M1->m, M1->n);
        mpz_mat_realloc(M2, d, d);
        mpz_mat_submat_swap(M2, 0, 0, M1, 0, 0, d, d);
        mpz_mat_determinant_triangular(det, M2);
        mpz_mod_ui(det, det, p);
        v.insert(std::make_pair(mpz_get_ui(det), M));
    }

    for(auto const & x : v) {
        printf("[det=%ld] ", x.first);
        mpz_mat_fprint(stdout, x.second);
    }

    mpz_clear(det);
}
