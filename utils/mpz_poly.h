#ifndef MPZ_POLY_H_
#define MPZ_POLY_H_

// IWYU pragma: no_include "double_poly.h"
// (only the fwd-decl is needed)
#include <stdio.h>
#include <stdint.h>
#include <gmp.h>
#ifdef __cplusplus
#include <sstream>
#include <string>
#include "fmt/core.h"
#endif
#include "macros.h"

#define xxxMPZ_POLY_TIMINGS
// for timings of roots mod p (beware, this is not thread-safe)

#ifndef DOUBLE_POLY_H_
typedef struct double_poly_s * double_poly_ptr;
typedef const struct double_poly_s * double_poly_srcptr;
#endif

#ifdef __cplusplus
#include <string>
#include <istream>      // std::istream // IWYU pragma: keep
#include <ostream>      // std::ostream // IWYU pragma: keep
extern "C" {
#endif

/* maximum degree we can reconstruct using mpz_poly_mul_tc_interpolate */
#define MAX_TC_DEGREE 19


/* Note, deg = -1 means P=0; otherwise, one should have coeff[deg] != 0.
   Warning: a polynomial of degree d needs d+1 allocation. */

struct mpz_poly_s {
  unsigned int alloc;
  int deg;
  mpz_t *coeff;
};
#ifndef DOUBLE_POLY_H_
/* double_poly.h forward-declares these. Don't do it twice */
typedef struct mpz_poly_s * mpz_poly_ptr;
typedef const struct mpz_poly_s * mpz_poly_srcptr;
#endif
typedef struct mpz_poly_s mpz_poly[1];

/* Note on parallelism.
 *
 * Some functions in this API can optionally use openmp. We want to make
 * sure that the code that calls these openmp-enabled functions does so
 * **willingly**. And we don't want to pollute the "normal" interface.
 *
 * The chosen solution is as follows.
 *
 * - the functions declared here, and use as plain (e.g.) mpz_poly_mul
 *   resolve to something that does **NOT** use openmp.
 *
 * - to use openmp, include mpz_poly_parallel.hpp instead, and call the
 *   member functions of an mpz_poly_parallel_info object. E.g.
 *
 *   mpz_poly_parallel_info inf;
 *   // (maybe add some configuration code for the inf object, if the
 *   // need for that ever appears)
 *   inf.mpz_poly_mul(....)
 *
 * More detail on can be found in mpz_poly_parallel.hpp and mpz_poly.cpp
 */

/* Management of the structure, set and print coefficients. */
void mpz_poly_init(mpz_poly_ptr, int d);
void mpz_poly_realloc (mpz_poly_ptr f, unsigned int nc);
void mpz_poly_set(mpz_poly_ptr g, mpz_poly_srcptr f);
void mpz_poly_swap (mpz_poly_ptr f, mpz_poly_ptr g);
void mpz_poly_clear(mpz_poly_ptr f);
static inline int mpz_poly_degree(mpz_poly_srcptr f) { return f->deg; }

void mpz_poly_cleandeg(mpz_poly_ptr f, int deg);
void mpz_poly_setcoeffs(mpz_poly_ptr f, mpz_t * coeffs, int d);
void mpz_poly_set_zero(mpz_poly_ptr f);
void mpz_poly_set_xi(mpz_poly_ptr f, int i);
void mpz_poly_set_mpz(mpz_poly_ptr f, mpz_srcptr z);
void mpz_poly_set_double_poly(mpz_poly_ptr g, double_poly_srcptr f);
/* returns 1 if parsing was successful */
int mpz_poly_set_from_expression(mpz_poly_ptr f, const char * value);

void mpz_poly_init_set_ab (mpz_poly_ptr rel, int64_t a, uint64_t b);
void mpz_poly_set_ab (mpz_poly_ptr rel, int64_t a, uint64_t b);
void mpz_poly_init_set_mpz_ab (mpz_poly_ptr rel, mpz_srcptr a, mpz_srcptr b);

void mpz_poly_setcoeff(mpz_poly_ptr f, int i, mpz_srcptr z);
void mpz_poly_setcoeff_si(mpz_poly_ptr f, int i, long z);
void mpz_poly_setcoeff_ui(mpz_poly_ptr f, int i, unsigned long z);
void mpz_poly_setcoeff_int64(mpz_poly_ptr f, int i, int64_t z);
void mpz_poly_setcoeff_uint64(mpz_poly_ptr f, int i, uint64_t z);
void mpz_poly_setcoeff_double(mpz_poly_ptr f, int i, double z);
void mpz_poly_getcoeff(mpz_ptr res, int i, mpz_poly_srcptr f);

/* functions for Joux-Lercier and generalized Joux-Lercier */
int mpz_poly_setcoeffs_counter(mpz_poly_ptr f, int* max_abs_coeffs, unsigned long *next_counter, int deg, unsigned long counter, unsigned int bound);
void  mpz_poly_setcoeffs_counter_print_error_code(int error_code);
unsigned long mpz_poly_getcounter(mpz_poly_ptr f, unsigned int bound);
unsigned long mpz_poly_cardinality(int deg, unsigned int bound);

/* return the leading coefficient of f */
static inline mpz_srcptr mpz_poly_lc (mpz_poly_srcptr f) {
    ASSERT(f->deg >= 0);
    return f->coeff[f->deg];
}

/* Print functions */
int mpz_poly_asprintf(char ** res, mpz_poly_srcptr f);
/* Print coefficients of f.
 * endl = 1 if "\n" at the end of fprintf. */
void mpz_poly_fprintf_endl (FILE *fp, mpz_poly_srcptr f, int endl);
void mpz_poly_fprintf(FILE *fp, mpz_poly_srcptr f);
void mpz_poly_fprintf_coeffs (FILE *fp, mpz_poly_srcptr f, const char sep);
void mpz_poly_fscanf_coeffs (FILE *fp, mpz_poly_ptr f, const char sep);
void mpz_poly_fprintf_cado_format (FILE *fp, mpz_poly_srcptr f,
                                   const char letter, const char *pre);
void mpz_poly_asprintf_cado_format (char **pstr, mpz_poly_srcptr f, const char letter,
                              const char *prefix);
void mpz_poly_print_raw(mpz_poly_srcptr f);
#ifdef MPZ_POLY_TIMINGS
  void print_timings_pow_mod_f_mod_p();
#endif
/* Tests and comparison functions */
int mpz_poly_cmp (mpz_poly_srcptr, mpz_poly_srcptr);
int mpz_poly_normalized_p (mpz_poly_srcptr f);
int mpz_poly_is_monic (mpz_poly_srcptr f);

/* Polynomial arithmetic */
void mpz_poly_to_monic(mpz_poly_ptr g, mpz_poly_srcptr f);
void mpz_poly_neg(mpz_poly_ptr f, mpz_poly_srcptr g);
void mpz_poly_add(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h);
void mpz_poly_sub(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h);
void mpz_poly_add_ui(mpz_poly_ptr g, mpz_poly_srcptr f, unsigned long a);
void mpz_poly_sub_ui(mpz_poly_ptr g, mpz_poly_srcptr f, unsigned long a);
void mpz_poly_add_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr a);
void mpz_poly_sub_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr a);
void mpz_poly_sub_mod_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h, mpz_srcptr m);
void mpz_poly_mul(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h);
void mpz_poly_mul_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_srcptr a);
void mpz_poly_divexact_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_srcptr a);
int mpz_poly_divisible_mpz (mpz_poly_srcptr P, mpz_srcptr a);
void mpz_poly_translation (mpz_poly_ptr, mpz_poly_srcptr, mpz_srcptr);
void mpz_poly_rotation (mpz_poly_ptr, mpz_poly_srcptr, mpz_poly_srcptr, mpz_srcptr, int);
void mpz_poly_addmul_si (mpz_poly_ptr, mpz_poly_srcptr, long);
void mpz_poly_mul_si (mpz_poly_ptr, mpz_poly_srcptr, long);
void mpz_poly_divexact_ui (mpz_poly_ptr, mpz_poly_srcptr, unsigned long);
void mpz_poly_rotation_int64 (mpz_poly_ptr, mpz_poly_srcptr, mpz_poly_srcptr, const int64_t, int);
void mpz_poly_makemonic_mod_mpz (mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_srcptr m);
void barrett_precompute_inverse (mpz_ptr invm, mpz_srcptr m);
int mpz_poly_mod_f_mod_mpz(mpz_poly_ptr R, mpz_poly_srcptr f, mpz_srcptr m, mpz_srcptr invf, mpz_srcptr invm);
int mpz_poly_mod_mpz(mpz_poly_ptr R, mpz_poly_srcptr A, mpz_srcptr m, mpz_srcptr invm);
int mpz_poly_mod_mpz_lazy (mpz_poly_ptr R, mpz_poly_srcptr A, mpz_srcptr m);
void mpz_poly_mul_mod_f_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P1, mpz_poly_srcptr P2, mpz_poly_srcptr f, mpz_srcptr m, mpz_srcptr invf, mpz_srcptr invm);
void mpz_poly_mul_mod_f (mpz_poly_ptr Q, mpz_poly_srcptr P1, mpz_poly_srcptr P2, mpz_poly_srcptr f);
void mpz_poly_reduce_frac_mod_f_mod_mpz (mpz_poly_ptr num, mpz_poly_ptr denom, mpz_poly_srcptr F, mpz_srcptr m);
int mpz_poly_div_qr (mpz_poly_ptr q, mpz_poly_ptr r, mpz_poly_srcptr f, mpz_poly_srcptr g, mpz_srcptr p);
int mpz_poly_div_r (mpz_poly_ptr h, mpz_poly_srcptr f, mpz_srcptr p);
int mpz_poly_div_qr_z (mpz_poly_ptr q, mpz_poly_ptr r, mpz_poly_srcptr f, mpz_poly_srcptr g);
int mpz_poly_div_r_z (mpz_poly_ptr r, mpz_poly_srcptr f, mpz_poly_srcptr g);
int mpz_poly_divexact (mpz_poly_ptr q, mpz_poly_srcptr h, mpz_poly_srcptr f, mpz_srcptr p);
void mpz_poly_div_2_mod_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr m);
void mpz_poly_div_xi(mpz_poly_ptr g, mpz_poly_srcptr f, int i);
void mpz_poly_mul_xi(mpz_poly_ptr g, mpz_poly_srcptr f, int i);
void mpz_poly_mul_xplusa(mpz_poly_ptr g, mpz_poly_srcptr f, mpz_srcptr a);

  
void mpz_poly_eval(mpz_ptr res, mpz_poly_srcptr f, mpz_srcptr x);
void mpz_poly_eval_ui (mpz_ptr res, mpz_poly_srcptr f, unsigned long x);
void mpz_poly_eval_diff_ui (mpz_ptr res, mpz_poly_srcptr f, unsigned long x);
void mpz_poly_eval_diff (mpz_ptr res, mpz_poly_srcptr f, mpz_srcptr x);
void mpz_poly_eval_poly(mpz_poly_ptr res, mpz_poly_srcptr f, mpz_poly_srcptr x);
void mpz_poly_eval_diff_poly (mpz_poly_ptr res, mpz_poly_srcptr f, mpz_poly_srcptr x);
void mpz_poly_eval_mod_mpz(mpz_t res, mpz_poly_srcptr f, mpz_srcptr x, mpz_srcptr m);
int mpz_poly_is_root(mpz_poly_srcptr poly, mpz_srcptr root, mpz_srcptr modulus);
void mpz_poly_eval_several_mod_mpz(mpz_ptr * res, mpz_poly_srcptr * f, int k, mpz_srcptr x, mpz_srcptr m);
void mpz_poly_sqr_mod_f_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, mpz_srcptr m, mpz_srcptr invf, mpz_srcptr invm);
void mpz_poly_pow_ui(mpz_poly_ptr B, mpz_poly_srcptr A, unsigned long n);
void mpz_poly_pow_ui_mod_f(mpz_poly_ptr B, mpz_poly_srcptr A, unsigned long n, mpz_poly_srcptr f);
void mpz_poly_pow_mod_f_mod_ui(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, mpz_srcptr a, unsigned long p);
void mpz_poly_pow_mod_f_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, mpz_srcptr a, mpz_srcptr p);
void mpz_poly_pow_ui_mod_f_mod_mpz (mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, unsigned long a, mpz_srcptr p);
void mpz_poly_derivative(mpz_poly_ptr df, mpz_poly_srcptr f);
mpz_poly* mpz_poly_base_modp_init (mpz_poly_srcptr P0, int p, unsigned long *K, int l);
void mpz_poly_base_modp_clear (mpz_poly *P, int l);
void mpz_poly_base_modp_lift(mpz_poly_ptr a, mpz_poly *P, int k, mpz_srcptr pk);
size_t mpz_poly_sizeinbase (mpz_poly_srcptr f, int base);
size_t mpz_poly_size (mpz_poly_srcptr f);
void mpz_poly_infinity_norm(mpz_ptr in, mpz_poly_srcptr f);
size_t mpz_poly_totalsize (mpz_poly_srcptr f);
void mpz_poly_gcd_mpz (mpz_poly_ptr h, mpz_poly_srcptr f, mpz_poly_srcptr g, mpz_srcptr p);
// compute f = GCD(f,g) mod N. If this fails, put the factor in the last
// given argument.
int mpz_poly_pseudogcd_mpz(mpz_poly_ptr , mpz_poly_ptr , mpz_srcptr , mpz_ptr);
void mpz_poly_xgcd_mpz(mpz_poly_ptr gcd, mpz_poly_srcptr f, mpz_poly_srcptr g, mpz_poly_ptr u, mpz_poly_ptr v, mpz_srcptr p);
void mpz_poly_homography (mpz_poly_ptr Fij, mpz_poly_srcptr F, int64_t H[4]);
void mpz_poly_homogeneous_eval_siui (mpz_ptr v, mpz_poly_srcptr f, const int64_t i, const uint64_t j);
void mpz_poly_content (mpz_ptr c, mpz_poly_srcptr F);
int mpz_poly_has_trivial_content (mpz_poly_srcptr F);
void mpz_poly_resultant(mpz_ptr res, mpz_poly_srcptr p, mpz_poly_srcptr q);
void mpz_poly_discriminant(mpz_ptr res, mpz_poly_srcptr f);
int mpz_poly_squarefree_p(mpz_poly_srcptr f);
int mpz_poly_is_irreducible_z(mpz_poly_srcptr f);

int mpz_poly_number_of_real_roots(mpz_poly_srcptr f);

struct mpz_poly_with_m_s {
    mpz_poly f;
    int m;
};
typedef struct mpz_poly_with_m_s mpz_poly_with_m[1];
typedef struct mpz_poly_with_m_s * mpz_poly_with_m_ptr;
typedef const struct mpz_poly_with_m_s * mpz_poly_with_m_srcptr;

struct mpz_poly_factor_list_s {
    mpz_poly_with_m * factors;
    int alloc;
    int size;
};
typedef struct mpz_poly_factor_list_s mpz_poly_factor_list[1];
typedef struct mpz_poly_factor_list_s * mpz_poly_factor_list_ptr;
typedef const struct mpz_poly_factor_list_s * mpz_poly_factor_list_srcptr;

void mpz_poly_factor_list_init(mpz_poly_factor_list_ptr l);
void mpz_poly_factor_list_clear(mpz_poly_factor_list_ptr l);
void mpz_poly_factor_list_flush(mpz_poly_factor_list_ptr l);
void mpz_poly_factor_list_push(mpz_poly_factor_list_ptr l, mpz_poly_srcptr f, int m);
void mpz_poly_factor_list_fprintf(FILE* fp, mpz_poly_factor_list_srcptr l);
int mpz_poly_factor_sqf(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f, mpz_srcptr p);
int mpz_poly_factor_ddf(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f0, mpz_srcptr p);
int mpz_poly_factor_edf(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f, int k, mpz_srcptr p, gmp_randstate_t rstate);

/* output is sorted by degree and lexicographically */
int mpz_poly_factor(mpz_poly_factor_list lf, mpz_poly_srcptr f, mpz_srcptr p, gmp_randstate_t rstate);
int mpz_poly_is_irreducible(mpz_poly_srcptr f, mpz_srcptr p);

/* lift from a factor list mod ell to a factor list mod ell2.
 * ell does not need to be prime, provided all factors considered are
 * unitary.
 *
 * ell and ell2 must be powers of the same prime, with ell2 <= ell^2
 * (NOTE that this is not checked)
 */
int mpz_poly_factor_list_lift(mpz_poly_factor_list_ptr fac, mpz_poly_srcptr f, mpz_srcptr ell, mpz_srcptr ell2);

/* This computes the ell-adic lifts of the factors of f, assuming
 * we have no multiplicities, using Newton lifting.
 * This requires that f be monic 
 *
 * I'm terribly lazy, so at the moment this is working only for prec==2.
 * Extending to arbitrary p is an easy exercise.
 *
 * The output is sorted based on the order of the factors mod p (that is,
 * factors are the lifts of the factors returned by mpz_poly_factor mod
 * p, in the same order).
 */
int mpz_poly_factor_and_lift_padically(mpz_poly_factor_list_ptr fac, mpz_poly_srcptr f, mpz_srcptr ell, int prec, gmp_randstate_t rstate);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/* This is sort of a generic way to write a c++ equivalent to the C type.
 * The first-class citizen in the cado-nfs code is (still) the C type, so
 * we're definitely bound to have a few infelicities here:
 *  - the type name can't be the same because of the size-1 array trick
 *    in C.
 *  - the C type is embedded as a member x for the same reason.
 *  - most operations on the C type should go through the member x
 *    (however, the conversions we have to _ptr and _srcptr can ease
 *    things a bit).
 */
struct cxx_mpz_poly {
    mpz_poly x;
    cxx_mpz_poly() { mpz_poly_init(x, -1); }
    inline int degree() const { return x->deg; } /* handy */
    cxx_mpz_poly(mpz_poly_srcptr f) { mpz_poly_init(x, -1); mpz_poly_set(x, f); }
    ~cxx_mpz_poly() { mpz_poly_clear(x); }
    cxx_mpz_poly(cxx_mpz_poly const & o) {
        mpz_poly_init(x, -1);
        mpz_poly_set(x, o.x);
    }
    cxx_mpz_poly & operator=(cxx_mpz_poly const & o) {
        mpz_poly_set(x, o.x);
        return *this;
    }
#if __cplusplus >= 201103L
    cxx_mpz_poly(cxx_mpz_poly && o) {
        mpz_poly_init(x, -1);
        mpz_poly_swap(x, o.x);
    }
    cxx_mpz_poly& operator=(cxx_mpz_poly && o) {
        mpz_poly_swap(x, o.x);
        return *this;
    }
#endif
    operator mpz_poly_ptr() { return x; }
    operator mpz_poly_srcptr() const { return x; }
    mpz_poly_ptr operator->() { return x; }
    mpz_poly_srcptr operator->() const { return x; }
    std::string print_poly(std::string const& var) const;
};

std::ostream& operator<<(std::ostream& o, cxx_mpz_poly const & f);
std::istream& operator>>(std::istream& in, cxx_mpz_poly & f);

#if GNUC_VERSION_ATLEAST(4,3,0)
extern void mpz_poly_init(cxx_mpz_poly & pl, int) __attribute__((error("mpz_poly_init must not be called on a mpz_poly reference -- it is the caller's business (via a ctor)")));
extern void mpz_poly_clear(cxx_mpz_poly & pl) __attribute__((error("mpz_poly_clear must not be called on a mpz_poly reference -- it is the caller's business (via a dtor)")));
#endif

struct mpz_poly_coeff_list {
    cxx_mpz_poly const & P;
    std::string sep;
    mpz_poly_coeff_list(cxx_mpz_poly const & P, std::string const & sep = ", "): P(P), sep(sep) {}
};
std::ostream& operator<<(std::ostream& os, mpz_poly_coeff_list const & P);

namespace fmt {
    template <> struct /* fmt:: */ formatter<mpz_poly_coeff_list>: formatter<string_view> {
    template <typename FormatContext>
        auto format(mpz_poly_coeff_list const & c, FormatContext& ctx) -> decltype(ctx.out())
        {
            std::ostringstream os;
            os << c;
            return formatter<string_view>::format( string_view(os.str()), ctx);
        }
};
    template <> struct /* fmt:: */ formatter<cxx_mpz_poly>: formatter<string_view> {
    template <typename FormatContext>
        auto format(cxx_mpz_poly const & c, FormatContext& ctx) -> decltype(ctx.out())
        {
            std::ostringstream os;
            os << c;
            return formatter<string_view>::format( string_view(os.str()), ctx);
        }
};
}


#endif

#endif	/* MPZ_POLY_H_ */
