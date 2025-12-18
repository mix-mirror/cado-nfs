#ifndef CADO_SIQS_SMALLSIEVE_GLUE_HPP
#define CADO_SIQS_SMALLSIEVE_GLUE_HPP

#include <cstdio>
#include <climits>

#ifdef HAVE_SSE41
#include <x86intrin.h>
#endif

#include "las-forwardtypes.hpp"         // spos_t
#include "fb-types.hpp"
#include "las-smallsieve-lowlevel.hpp"

#include "sieve-methods.hpp"

#include "macros.h"

struct list_nil {};

template<typename T, int b, typename F> struct choice_list_car {};

struct siqs_small_sieve_best_code_choices {
#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)

    /* Mostly copied from las-smallsieve-glue.hpp.
     * TODO needs to be benchmarked
     */
    typedef choice_list_car<
                assembly_generic_oldloop,
                5, /* use the function above for primes that have a max number
                      of hits that is at least 32 */
            choice_list_car<
                assembly4,
                4, /* use the function above for primes that have a max number
                      of hits that is at least 16 */
            choice_list_car<
                assembly3,
                3, /* use the function above for primes that have a max number
                      of hits that is at least 8 */
            choice_list_car<
                assembly2x,
                0, /* use the function above for the remaining primes */
            list_nil
            >>>> type;
#else
    typedef choice_list_car<
                manual_oldloop,
                0,
            list_nil
            > type;
#endif
};

/* small sieve classes */
/* So many things are used in common for many small sieve routines that
 * it makes sense to gather them in a common object */
struct siqs_small_sieve_base {
    int min_logI_logB;
    int logI;
    int N;
    unsigned int log_regions_per_line;
    unsigned int region_rank_in_line;
    bool last_region_in_line;
    unsigned int j0;
    unsigned int j1;
    int i0;
    inline int F() const { return 1 << min_logI_logB; }
    inline int I() const { return 1 << logI; }
    siqs_small_sieve_base(int logI, int N)
        : logI(logI), N(N)
    {
        unsigned int log_lines_per_region;
        min_logI_logB         = std::min(LOG_BUCKET_REGION, logI);
        log_lines_per_region  = (LOG_BUCKET_REGION-min_logI_logB);
        log_regions_per_line  = (logI-min_logI_logB);
        unsigned int regions_per_line      = (1<<log_regions_per_line);
        region_rank_in_line   = (N&(regions_per_line-1));
        j0    = ((N>>log_regions_per_line)<<log_lines_per_region);
        j1    = (j0+(1<<log_lines_per_region));
        i0    = ((region_rank_in_line<<LOG_BUCKET_REGION)-(1 << (logI-1)));
    }

    spos_t first_position_ordinary_prime(
            siqs_ssp_simple_t const & ssp,
            unsigned int dj) const
    {
        fbprime_t p = ssp.get_p();
        fbprime_t r = ssp.get_r(); /* it corresponds to the root for j = 0 */
        unsigned int j = j0 + dj;
        for (auto const rk: ssp.CRT_data()) {
            if (j & 1u) {
                r = addmod_u32(r, rk, p);
            }
            j >>= 1u;
        }
        int64_t x = (int64_t) r - i0;
        x = x % (int64_t) p;
        if (x < 0) x += p;
        return (spos_t) x;
    }

    spos_t first_position_power_of_two(
            siqs_ssp_t const & ssp,
            unsigned int dj) const
    {
        fbprime_t pmask = ssp.get_p() - 1u;
        fbprime_t r = ssp.get_r(); /* it corresponds to the root for j = 0 */
        unsigned int j = j0 + dj;
        for (auto const rk: ssp.CRT_data()) {
            if (j & 1u) {
                r = (r + rk) & pmask;
            }
            j >>= 1u;
        }
        int64_t x = (int64_t) r - i0;
        x = x & pmask;
        if (x < 0) x += ssp.get_p();
        return (spos_t) x;
    }

    fbroot_t first_position_in_line(
            siqs_ssp_t const & ssp,
            unsigned int dj) const
    {
        fbprime_t p = ssp.get_p();
        fbprime_t r = ssp.get_r(); /* it corresponds to the root for j = 0 */
        unsigned int j = j0 + dj;
        for (auto const rk: ssp.CRT_data()) {
            if (j & 1u) {
                r = addmod_u32(r, rk, p);
            }
            j >>= 1u;
        }
        int64_t x = (int64_t) r - i0;
        x = x % (int64_t) p;
        if (x < 0) x += p;
        return (spos_t) x;
    }

};

struct siqs_small_sieve : public siqs_small_sieve_base {
    typedef siqs_small_sieve_base super;
    std::vector<int> const & positions;
    std::vector<siqs_ssp_simple_t> const& primes;
    size_t sorted_limit;
    std::list<size_t> sorted_subranges;
    std::vector<siqs_ssp_t> const& not_nice_primes;
    unsigned char * S;

    static const fbprime_t pattern2_size = 2 * sizeof(unsigned long);
    static const int test_divisibility = 0; /* very slow, but nice for debugging */

    /* This "index" field is actually the loop counter that is shared by
     * the various instantiations that are triggered by
     * handle_nice_primes_meta_loop<>.
     */
    size_t index = 0;

    siqs_small_sieve(
            std::vector<int> const& positions,
            std::vector<siqs_ssp_simple_t> const & primes,
            std::vector<siqs_ssp_t> const & not_nice_primes,
            unsigned char*S, int logI,
            unsigned int N)
        : super(logI, N)
        , positions(positions)
        , primes(primes)
        , not_nice_primes(not_nice_primes)
        , S(S)
    {
        auto s = begin(primes);
        auto s0 = s;
        auto e = end(primes);
        for( ; e > s ; ) {
            if ((e-s) < 32) {
                /* we don't want to bother adding an extra
                 * control loop for a small bunch of primes
                 */
                break;
            }
            int c;
            for(c = 1 ; c < (e-s) && !(s[c] < s[c-1]) ; c++);
            if (c <= 16) {
                if (s != s0) {
                    /* s != s0 : cf mail ET->{PG,PZ} 201906071825
                     */
                    fprintf(stderr, "warning, the prime list looks really ugly\n");
                }
            }
            /*
            fprintf(stderr, "ssp entries [%zd..%zd[ (out of %zu) are sorted\n",
                    s-s0, s+c-s0, primes.size());
                    */
            sorted_subranges.push_back((s+c) - s0);
            s += c;
        }
        /*
        fprintf(stderr, "ssp entries [%zd..%zd[ (tail) do not need to be sorted\n",
                s-s0, primes.size());
                */
    }

    bool finished() const { return index == primes.size(); }
    bool finished_sorted_prefix() const { return index == sorted_limit; }

    /* most fields from the parent class must be redeclared here if we
     * want to use them transparently -- it's only a shorthand
     * convenience. We don't re-expose everything though, since some are
     * only seldom used.
     */
    using super::logI;
    using super::i0;
    using super::j0;
    using super::j1;
    using super::I;
    using super::F;

    void handle_power_of_2(
            siqs_ssp_t const & ssp,
            where_am_I & w MAYBE_UNUSED)
    {
        /* Powers of 2 are treated separately */
        /* Don't sieve powers of 2 again that were pattern-sieved */
        const fbprime_t p = ssp.get_p();
        WHERE_AM_I_UPDATE(w, p, p);

        if (ssp.is_pattern_sieved())
            return;

        const unsigned char logp = ssp.logp;
        unsigned char *S_ptr = S;

        for (unsigned int dj = 0; dj < j1-j0; ++dj) {
            //TODO be more efficient (gray code)
            int pos = super::first_position_power_of_two(ssp, dj);
            for (int i = pos; i < F(); i += p) {
                WHERE_AM_I_UPDATE(w, x, ((size_t) dj << logI) + i);
                sieve_increase(S_ptr + i, logp, w);
            }
            S_ptr += I();
        }
    }

    template<typename inner_loop, int bits_off>
    bool handle_nice_prime(
            siqs_ssp_simple_t const & ssp,
            MAYBE_UNUSED spos_t pos,
            where_am_I & w)
    {
        const fbprime_t p = ssp.get_p();
        if (bits_off != 0 && (p >> (super::min_logI_logB + 1 - bits_off))) {
            /* time to move on to the next bit size; */
            return false;
        }

        MAYBE_UNUSED const fbprime_t r = ssp.get_r();
        WHERE_AM_I_UPDATE(w, r, r);
        const unsigned char logp = ssp.logp;
        unsigned char * S0 = S;
        unsigned char * S1 = S + F();

        /* we sieve over the area [S0..S0+(i1-i0)] (F() is (i1-i0)),
         * which may actually be just a fragment of a line. After
         * that, if (i1-i0) is different from I, we'll break anyway.
         * So whether we add I or (i1-i0) to S0 does not matter much.
         */

        //TODO be more efficient (gray code)
        for (unsigned int dj = 0, j = j0; j < j1; ++dj, ++j) {
            WHERE_AM_I_UPDATE(w, j, dj);
            int pos = super::first_position_ordinary_prime(ssp, dj);
            inner_loop()(S0, S1, S0 - S, pos, p, logp, w);
            S0 += I();
            S1 += I();
        }
        return true;
    }

    /* This function is responsible for small-sieving primes of a
     * specific bit size */
    template<typename inner_loop, int bits_off>
    inline void handle_nice_primes(where_am_I & w MAYBE_UNUSED)
    {
        /* here, we can sieve for primes p < 2 * F() / 2^bits_off,
         * (where F() is i1-i0 = 2^min(logI, logB)).
         *
         * meaning that the number of hits in a line is at least
         * floor(F() / p) = 2^(bits_off-1)
         *
         * Furthermore, if p >= 2 * F() / 2^(bits_off+1), we can also
         * say that the number of hits is at most 2^bits_off
         */

        for( ; index < sorted_limit ; index++) {
            auto const & ssp(primes[index]);
            spos_t pos = positions[index];
            const fbprime_t p = ssp.get_p();
            if (bits_off != 0 && (p >> (super::min_logI_logB + 1 - bits_off))) {
                /* time to move on to the next bit size; */
                return;
            }
            WHERE_AM_I_UPDATE(w, p, p);
            handle_nice_prime<inner_loop, bits_off>(ssp, pos, w);
        }
    }

    private:
    /* template machinery to make one single function out of several */
    /* we'll now craft all the specific handle_nice_primes_meta_loop
     * functions into one big function. Because we're playing tricks with
     * types and lists of types and such, we need to work with partial
     * specializations at the class level, which is admittedly messy. */
    template<typename T, int max_bits_off = INT_MAX>
        struct handle_nice_primes_meta_loop
        {
            void operator()(siqs_small_sieve & SS, where_am_I &) {
                /* default, should be at end of list. We require that we
                 * are done processing, at this point. */
                ASSERT_ALWAYS(SS.finished_sorted_prefix());
            }
        };

    /* optimization: do not split into pieces when we have several times
     * the same code anyway. */
    template<typename L0, int b0, int b1, typename T, int bn>
        struct handle_nice_primes_meta_loop<choice_list_car<L0,b0,
                     choice_list_car<L0,b1,
                    T>>, bn>
        {
            static_assert(b0 > b1, "choice list is in wrong order");
            void operator()(siqs_small_sieve & SS, where_am_I & w) {
                handle_nice_primes_meta_loop<choice_list_car<L0,b1, T>, bn>()(SS, w);
            }
        };
    template<typename L0, int b0, int bn>
        struct is_compatible_for_range {
            static_assert(bn > b0, "choice list is in wrong order");
            static const int value =
                L0::template is_compatible<bn>::value &&
                is_compatible_for_range<L0, b0, bn-1>::value;
        };
    template<typename L0, int b0>
        struct is_compatible_for_range<L0, b0, b0> {
            static const int value =
                L0::template is_compatible<b0>::value;
        };
    template<typename L0, int b0>
        struct is_compatible_for_range<L0, b0, INT_MAX> {
            static const int value =
                L0::template is_compatible<INT_MAX>::value &&
                is_compatible_for_range<L0, b0, b0 + 5>::value;
        };

    template<typename L0, int b0, typename T, int bn>
        struct handle_nice_primes_meta_loop<choice_list_car<L0,b0,T>, bn>
        {
            static_assert(is_compatible_for_range<L0, b0, bn>::value, "Cannot use this code fragment for primes of the current size");
            void operator()(siqs_small_sieve & SS, where_am_I & w) {
                SS.handle_nice_primes<L0, b0>(w);
                handle_nice_primes_meta_loop<T, b0-1>()(SS, w);
            }
        };

    public:
    void do_pattern_sieve(where_am_I &);

    void normal_sieve(where_am_I & w) {
        for(size_t s : sorted_subranges) {
            sorted_limit = s;
            /* This function will eventually call handle_nice_primes on
             * sub-ranges of the set of small primes. */
            typedef siqs_small_sieve_best_code_choices::type choices;
            siqs_small_sieve::handle_nice_primes_meta_loop<choices>()(*this, w);
        }

        /* This is for the tail of the list. We typically have powers,
         * here. These are ordinary, nice, simple prime powers, but the
         * only catch is that these don't get resieved (because the prime
         * itself was already divided out, either via trial division or
         * earlier resieving. By handling them here, we benefit from the
         * ssdpos table. */
        for( ; index < primes.size() ; index++) {
            auto const & ssp(primes[index]);
            spos_t pos = positions[index];
            WHERE_AM_I_UPDATE(w, p, ssp.get_p());
            handle_nice_prime<default_smallsieve_inner_loop, 0>(ssp, pos, w);
        }
    }

    void exceptional_sieve(where_am_I & w) {
        /* a priori we'll never have "nice" primes here, but we're not
         * forced to rule it out completely, given that we have the code
         * available at hand. The only glitch is that we're storing the
         * start positions *only* for the primes in the ssps array. */

        for(auto const & ssp : not_nice_primes) {
            if (ssp.is_pattern_sieved()) {
                /* This ssp is pattern-sieved, nothing to do here */
            } else if (ssp.is_pow2()) {
                handle_power_of_2(ssp, w);
            } else {
                /* I don't think we can end up here.  */
                ASSERT_ALWAYS(0);
            }
        }
    }

};

#endif	/* CADO_SIQS_SMALLSIEVE_GLUE_HPP */
