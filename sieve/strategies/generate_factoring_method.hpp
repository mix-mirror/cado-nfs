#ifndef GENERATE_FACTORING_METHOD_HPP
#define GENERATE_FACTORING_METHOD_HPP

#include <istream>
#include <ostream>
#include <vector>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "facul.hpp"
#include "facul_ecm.h"
#include "tab_fm.hpp"
#include "tab_point.hpp"

/************************************************************************/
/*                      To generate numbers n=pq*/
/************************************************************************/

double *distribution_prime_number(int min, int max);

cxx_mpz generate_prime_factor(gmp_randstate_t state, int lenFact);

cxx_mpz
generate_composite_integer(gmp_randstate_t state,
			   int lenFact1, int lenFactall);

int
select_random_index_according_dist(double const * dist, int len,
                                   gmp_randstate_ptr state);

cxx_mpz
generate_composite_integer_interval(gmp_randstate_t state,
				    double const * dist, int lenFact1_min,
				    int lenFact1_max, int lenFactall);

/************************************************************************/
/*                To model our factoring methods                        */
/************************************************************************/

/* This generate a 1-sided strategy with just one single method, and (for
 * ECM) a parameter sigma drawn from state. extra_primes is set to 1.
 */
facul_strategy_oneside generate_fm (facul_method_code method,
                               unsigned long B1, unsigned long B2,
                               ec_parameterization_t curve,
                               gmp_randstate_ptr state
                               );

/************************************************************************/
/*                      Create Factoring Methods */
/************************************************************************/

int *choice_parameters(facul_method_code method, int len_p_min);

tabular_fm bench_proba_time_pset(facul_method_code method, ec_parameterization_t curve,
				    gmp_randstate_t state,
				    int len_p_min, int len_p_max,
				    int len_n, const int *param_region);

tabular_fm generate_factoring_methods(gmp_randstate_t state, int len_p_min,
					 int len_p_max, int len_n, int opt_ch,
					 const int *param_sieve);

tabular_fm generate_factoring_methods_mc(gmp_randstate_t state,
					    int len_p_min, int len_p_max,
					    int len_n, facul_method_code method, ec_parameterization_t curve,
					    int opt_ch,
                                            const int *param_sieve);

/************************************************************************/
/*                      ANALYSE AND FILTER */
/************************************************************************/

double bench_proba_fm(facul_strategy_oneside const & strategy,
        gmp_randstate_t state,
 		      unsigned long len_p, unsigned long len_n,
                      std::vector<cxx_mpz> & N,
		      size_t nb_test_max);

void bench_proba(gmp_randstate_t state, tabular_fm & fm, int len_p_min,
        int p_max, size_t nb_test_max);


double bench_time_fm_onelength(facul_strategy_oneside const & method,
                               std::vector<cxx_mpz> & N,
			       size_t nb_test);

void bench_time(gmp_randstate_t state, tabular_fm & fm, size_t nb_test);


tabular_fm filtering(tabular_fm const & fm, int final_nb_methods);

/************************************************************************/
/*                      CONVEX_HULL_FM                                  */
/************************************************************************/

tabular_fm convex_hull_fm(tabular_fm const & t);

/************************************************************************/
/*               COMPUTE THE CONVEX_HULL FROM A FILE OF FM              */
/************************************************************************/

tabular_fm convex_hull_from_file(std::istream & in, std::ostream & out);

struct weighted_success {
    double prob = 0;
    double time = 0;
    // weighted_success() = default;
    weighted_success(double p, double t) : prob(p), time(t) {}
    weighted_success(size_t p, double t, size_t N)
        : prob(double(p)/double(N))
        , time(t/double(N))
    {}
};

#endif				/* GENERATE_FACTORING_METHOD_HPP */
