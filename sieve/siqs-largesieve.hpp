#ifndef CADO_SIQS_LARGESIEVE_HPP
#define CADO_SIQS_LARGESIEVE_HPP

#include <cstdint>

#include <algorithm>
#include <limits>
#include <tuple>
#include <vector>

#include "las-arith.hpp"
#include "las-qlattice.hpp"
#include "macros.h"

class siqs_largesieve
{
  public:
    using T_elt = std::pair<uint32_t, uint32_t>;

  protected:
    uint32_t pp;
    uint32_t invq; /* 1/q modulo pp */
    /* crt_data_modp := [ Rk/q mod p for Rk in Q.crt_data_modq ] */
    std::vector<uint32_t> crt_data_modp;
    std::vector<uint32_t> roots; /* r/q modulo pp, for r roots modulo pp */
    size_t n; /* number of factors in the special-q */
    size_t n1, n2;
    std::vector<T_elt> T1, T2;

    void set_Ti(std::vector<T_elt> & T, size_t start, size_t end)
    {
        T = {{0U, 0U}};
        std::vector<T_elt> Tp, Tm;
        uint32_t mask = 1U << start;

        auto comp = [](T_elt const & u, T_elt const & v) -> bool {
            return u.first < v.first;
        };

        for (size_t k = start; k < end; ++k, mask<<=1U) {
            size_t new_size = 2*T.size();
            shift(Tp, T, crt_data_modp[k], mask);
            shift(Tm, T, pp-crt_data_modp[k], 0U);

            T.clear();
            T.reserve(new_size);
            std::ranges::merge(Tp, Tm, std::back_inserter(T), comp);
            Tp.clear();
            Tm.clear();
        }
    }

  public:
    template<class FB_ENTRY_TYPE>
    siqs_largesieve(
            siqs_special_q_data const & Q,
            FB_ENTRY_TYPE const & e,
            size_t n)
        : pp(e.get_q())
        , n(n)
        , n1(n/2)
        , n2(n-n1)
    {
        ASSERT_ALWAYS(n <= Q.nfactors());
        e.compute_crt_data_modp(invq, crt_data_modp, Q, false);
        /* copy roots */
        if constexpr (std::is_same_v<FB_ENTRY_TYPE, fb_entry_general>)
        {
            for(unsigned char i=0U; i < e.nr_roots; ++i) {
                roots.push_back(e.roots[i].r);
            }
        } else {
            for (auto const & r: e.roots) {
                roots.push_back(r);
            }
        }
        /* transform them: r/q mod pp */
        for (auto & r: roots) {
            if (UNLIKELY(!(pp & 1))) { /* p power of 2 */
                r = (r * invq) & (pp - 1U);
            } else {
                r = mulmodredc_u32<true>(r, invq, pp, e.invq);
            }
        }

        set_Ti(T1, 0U, n1);
        set_Ti(T2, n1, n1+n2);
    }

    /* Assumes 0 <= r < pp */
    void shift(std::vector<T_elt> & Tt, std::vector<T_elt> const & T,
               uint32_t const r, uint32_t const mask)
    {
        /* despite the name, it performs a binary search */
        auto zp = std::ranges::lower_bound(T, pp-r, {}, &T_elt::first);

        Tt.reserve(T.size());
        auto op1 = [r, mask, pp = this->pp](T_elt const & e) -> T_elt {
            return { e.first + r - pp, e.second | mask};
        };
        auto op2 = [r, mask](T_elt const & e) -> T_elt {
            return { e.first + r, e.second | mask};
        };
        std::ranges::transform(zp, T.end(), std::back_inserter(Tt), op1);
        std::ranges::transform(T.begin(), zp, std::back_inserter(Tt), op2);
    }

    void prepare_for_root(
            std::vector<T_elt> & T1s,
            size_t const root_idx,
            int const logI,
            uint32_t const highest_bits)
    {
        T1s.clear();
        uint32_t s = (roots[root_idx] + (1 << (logI-1))) % pp;
        uint32_t h = highest_bits;
        for (size_t i = n; i < crt_data_modp.size(); ++i, h>>=1) {
            if (h & 1u) {
                s = (s + crt_data_modp[i]) % pp;
            } else {
                s = (s + (pp - crt_data_modp[i])) % pp;
            }
        }

        shift(T1s, T1, s, highest_bits << n);
    }

    uint32_t get_pp() const {
        return pp;
    }

    std::vector<T_elt> const & get_T2() const {
        return T2;
    }
};

#endif /* CADO_SIQS_LARGESIEVE_HPP */
