#include "cado.h" // IWYU pragma: keep

#include <cfloat>
#include <cmath>

#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "fmt/base.h"

#include "convex_hull.hpp"
#include "decomp.hpp"
#include "facul_ecm.h"
#include "facul_method.hpp"
#include "fm.hpp"
#include "generate_strategies.hpp"
#include "macros.h"
#include "point.hpp"
#include "strategy.hpp"
#include "tab_decomp.hpp"
#include "tab_fm.hpp"
#include "tab_point.hpp"
#include "tab_strategy.hpp"
#include "utils_cxx.hpp"

/* Chaining PM1 after PP1 (or the other way round) does not give
 * independent events, so the failure probabilities do not simply
 * multiply. This is the correction that has always been applied. */
static bool is_pm1_or_pp1(factoring_method const & fm)
{
    return fm.params.method == PM1_METHOD ||
           fm.params.method == PP1_27_METHOD ||
           fm.params.method == PP1_65_METHOD;
}

/************************************************************************/
/*                      COLLECT DATA FOR ONLY ONE COFACTOR              */
/************************************************************************/

/* The probability that 'fm' fails to find a non-trivial factor of a
 * number with this decomposition. */
double compute_proba_method_one_decomp(decomp const & D,
                                       factoring_method const & fm)
{
    double proba_fail = 1;
    for (auto j: D) {
        j -= fm.len_p_min;
        if (j < fm.proba.size())
            proba_fail *= 1 - fm.proba[j];
        // else the probability is close enough to 0.
    }
    return proba_fail;
}

/* The probability that 'strat' finds a non-trivial factor, averaged over
 * the decompositions that fall within [len_p_min, len_p_max]. */
double compute_proba_strategy(tabular_decomp const & init_tab,
                              strategy_t const & strat,
                              unsigned int len_p_min, unsigned int len_p_max)
{
    double all = 0.0;
    double nb_found_elem = 0;

    for (auto const & D: init_tab) {
        if (is_good_decomp(D, len_p_min, len_p_max)) {
            double p_fail_all = 1;
            for (auto const & fm: strat.tab_fm) {
                double const p_fail_one = compute_proba_method_one_decomp(D, fm);
                p_fail_all *= p_fail_one;
                if (is_pm1_or_pp1(fm))
                    p_fail_all = (p_fail_one + p_fail_all) / 2;
            }
            nb_found_elem += (1 - p_fail_all) * D.nb_elem;
        }
        all += D.nb_elem;
    }
    if (all < (double)LDBL_EPSILON) // no decomposition at all
        return 0;
    return nb_found_elem / all;
}

/* The average time 'strat' takes on a cofactor of r bits. */
double compute_time_strategy(tabular_decomp const & init_tab,
                             strategy_t const & strat, unsigned int r)
{
    double time_average = 0;
    double all = 0.0;

    for (auto const & D: init_tab) {
        double time_dec = 0;
        double proba_fail_all = 1;
        for (auto const & fm: strat.tab_fm) {
            time_dec += fm.time_for(r) * proba_fail_all;

            double const proba_fail_method =
                compute_proba_method_one_decomp(D, fm);
            proba_fail_all *= proba_fail_method;
            if (is_pm1_or_pp1(fm))
                proba_fail_all = (proba_fail_all + proba_fail_method) / 2;
        }

        time_average += time_dec * D.nb_elem;
        all += D.nb_elem;
    }
    if (all < (double)LDBL_EPSILON) // no decomposition at all
        return 0;

    return time_average / all;
}

/* Append 'strategy' to 't', dropping the zero methods it holds. A
 * strategy made only of zero methods keeps its first one, so that it is
 * not empty. */
static void add_strategy_without_zero(tabular_strategy & t,
                                      strategy_t const & strategy)
{
    strategy_t elem;
    bool strat_is_zero = true;
    for (auto const & fm: strategy.tab_fm) {
        if (!fm.is_zero()) {
            strat_is_zero = false;
            elem.add_fm(fm);
        }
    }
    if (strat_is_zero)
        elem.add_fm(strategy.tab_fm[0]);

    elem.proba = strategy.proba;
    elem.time = strategy.time;
    t.push_back(std::move(elem));
}

/************************************************************************/
/*                   GENERATE MATRIX                                    */
/************************************************************************/

/* Add the different chains of ECM curves to 'strat'. 'lbucket' lets
 * curves be added in buckets of that length. */
static void generate_collect_iter_ecm(
    tabular_fm const & ecm, unsigned int ind_ecm, strategy_t & strat,
    unsigned int ind_tab, unsigned int index_iter, unsigned int len_iteration,
    unsigned int lbucket, tabular_decomp const & init_tab,
    tabular_strategy & all_strat, unsigned int fbb, unsigned int lpb,
    unsigned int r, int is_already_used_B12M16)
{
    if (index_iter >= len_iteration) {
        add_strategy_without_zero(all_strat, strat);
        strategy_t & added = all_strat.back();
        added.proba = compute_proba_strategy(init_tab, added, fbb, lpb);
        added.time = compute_time_strategy(init_tab, added, r);
        return;
    }

    for (unsigned int i = ind_ecm; i < ecm.size(); i++) {
        /* BRENT12 and MONTY16 have only one useful sigma, so they are
         * used at most once in a chain. */
        if (ecm[i].params.parameterization == MONTY16 ||
            ecm[i].params.parameterization == BRENT12) {
            if (is_already_used_B12M16)
                continue;
            strat.tab_fm[ind_tab] = ecm[i];
            generate_collect_iter_ecm(ecm, i + 1, strat, ind_tab + 1,
                                      index_iter + 1, len_iteration, lbucket,
                                      init_tab, all_strat, fbb, lpb, r, true);
        } else { // MONTY12
            for (unsigned int j = 0;
                 j < lbucket && ind_tab + j < strat.tab_fm.size(); j++)
                strat.tab_fm[ind_tab + j] = ecm[i];
            generate_collect_iter_ecm(ecm, i, strat, ind_tab + lbucket,
                                      index_iter + lbucket, len_iteration,
                                      lbucket + 1, init_tab, all_strat, fbb,
                                      lpb, r, is_already_used_B12M16);
        }
    }

    /* Keep the RAM in check by reducing the collection to its convex
     * hull. */
    if (all_strat.size() < 100000)
        all_strat = convex_hull_strategy(all_strat);
}

/* The best strategies for one cofactor size 'r' on one side. The
 * generator chains methods as
 *   PM1 (0/1) + PP1 (0/1) + ECM-M12(0/1/2...) + ECM-M16/B12(0/1) + ECM-M12(...)
 * and keeps only the convex hull as it goes, so that the collection does
 * not fill the RAM.
 */
tabular_strategy generate_strategies_oneside(tabular_decomp const & init_tab,
                                             factoring_method const & zero,
                                             tabular_fm const & pm1,
                                             tabular_fm const & pp1,
                                             tabular_fm const & ecm,
                                             int ncurves, unsigned long lim,
                                             unsigned int lpb, unsigned int r)
{
    unsigned int const fbb = ceil(log2((double)(lim + 1)));
    unsigned int const lim_is_prime = 2 * fbb - 1;

    ASSERT_ALWAYS(!init_tab.empty() == (r >= lim_is_prime));

    if (r < lim_is_prime) {
        /* r is already prime. Two zero strategies cover the cases: a
         * good prime (fbb < r < lpb) succeeds with probability 1, a
         * prime that is too large, or a length that cannot occur, with
         * probability 0. */
        strategy_t st_zero;
        st_zero.add_fm(zero);
        st_zero.time = 0.0;
        st_zero.proba = (r != 1 && (r < fbb || r > lpb)) ? 0 : 1.0;

        tabular_strategy res;
        res.push_back(std::move(st_zero));
        return res;
    }

    tabular_strategy all_strat;

    int const len_strat = 2 + ncurves;
    strategy_t strat;
    for (int i = 0; i < len_strat; i++)
        strat.add_fm(zero);

    for (auto const & fm_pm1: pm1) {
        double const current_proba_pm1 = fm_pm1.proba[0];
        strat.tab_fm[0] = fm_pm1;
        for (auto const & fm_pp1: pp1) {
            if (fm_pp1.params.B1 != 0 && fm_pp1.proba[0] < current_proba_pm1)
                continue;
            strat.tab_fm[1] = fm_pp1;

            generate_collect_iter_ecm(ecm, 0, strat, 2, 0, ncurves, 0, init_tab,
                                      all_strat, fbb, lpb, r, false);
        }
    }
    return convex_hull_strategy(all_strat);
}

/* Concatenate two one-sided strategies, st1 going first on
 * 'first_side'. */
static strategy_t concat_strategies(strategy_t const & st1,
                                    strategy_t const & st2, int first_side)
{
    strategy_t st;
    for (auto const & fm: st1.tab_fm)
        st.add_fm(fm, first_side);
    for (auto const & fm: st2.tab_fm)
        st.add_fm(fm, first_side ? 0 : 1);
    return st;
}

/* The best strategies for a pair of cofactors, given the optimal
 * strategies for each side. The probability and time of each side must
 * already have been computed. */
tabular_strategy generate_strategy_r0_r1(tabular_strategy const & strat_r0,
                                         tabular_strategy const & strat_r1)
{
    tabular_strategy strat_r0_r1;
    tabular_strategy ch;
    unsigned long nb_strat = 0;

    for (size_t r = 0; r < strat_r0.size(); r++) {
        double const p0 = strat_r0[r].proba;
        double const c0 = strat_r0[r].time;
        for (auto const & s1: strat_r1) {
            nb_strat++;
            double const p1 = s1.proba;
            double const c1 = s1.time;
            double const proba = p0 * p1;
            double const tps0 = c0 + p0 * c1; // starting with side 0
            double const tps1 = c1 + p1 * c0; // starting with side 1
            strategy_t st;
            if (tps0 < tps1) {
                st = concat_strategies(strat_r0[r], s1, 0);
                st.time = tps0;
            } else {
                st = concat_strategies(s1, strat_r0[r], 1);
                st.time = tps1;
            }
            st.proba = proba;
            strat_r0_r1.push_back(std::move(st));
        }

        /* Keep the RAM in check. */
        if (nb_strat > 100000 || r == strat_r0.size() - 1) {
            strat_r0_r1.insert(strat_r0_r1.end(), ch.begin(), ch.end());
            ch = convex_hull_strategy(strat_r0_r1);
            strat_r0_r1.clear();
            nb_strat = 0;
        }
    }

    return ch;
}

/* The best strategies for each pair of cofactor sizes, from a set of
 * factoring methods whose probabilities and times have already been
 * measured (with gfm). */
strategy_matrix generate_matrix(std::string const & name_directory_decomp,
                                tabular_fm const & pm1, tabular_fm const & pp1,
                                tabular_fm const & ecm, int ncurves,
                                unsigned long lim0, unsigned int lpb0,
                                unsigned int mfb0, unsigned long lim1,
                                unsigned int lpb1, unsigned int mfb1)
{
    unsigned int const fbb0 = ceil(log2((double)(lim0 + 1)));
    unsigned int const fbb1 = ceil(log2((double)(lim1 + 1)));

    strategy_matrix matrix(mfb0 + 1);
    for (auto & row: matrix)
        row.resize(mfb1 + 1);

    factoring_method const zero = factoring_method::from_fields(0, 0, 0, 0);

    auto read_decomp = [&](unsigned long lim, unsigned int r,
                           unsigned int lim_is_prime) {
        tabular_decomp tab_decomp;
        if (r >= lim_is_prime) {
            auto filename =
                fmt::format("{}/decomp_{}_{}", name_directory_decomp, lim, r);
            std::ifstream is(filename);
            if (!(is >> tab_decomp))
                throw cado::error("Cannot read {}", filename);
        }
        return tab_decomp;
    };

    /* Side 0 is precomputed for every size; side 1 is computed as we go. */
    std::vector<tabular_strategy> data_rat(mfb0 + 1);
    for (unsigned int r0 = 0; r0 <= mfb0; r0++) {
        auto const tab_decomp = read_decomp(lim0, r0, 2 * fbb0 - 1);
        data_rat[r0] = generate_strategies_oneside(tab_decomp, zero, pm1, pp1,
                                                   ecm, ncurves, lim0, lpb0, r0);
    }

    for (unsigned int r1 = 0; r1 <= mfb1; r1++) {
        auto const tab_decomp = read_decomp(lim1, r1, 2 * fbb1 - 1);
        tabular_strategy const strat_r1 = generate_strategies_oneside(
            tab_decomp, zero, pm1, pp1, ecm, ncurves, lim1, lpb1, r1);

        for (unsigned int r0 = 0; r0 <= mfb0; r0++)
            matrix[r0][r1] = generate_strategy_r0_r1(data_rat[r0], strat_r1);
    }

    return matrix;
}

/************************************************************************/
/*                      CONVEX_HULL_ST                                  */
/************************************************************************/

tabular_point convert_tab_point_to_tab_strategy(tabular_strategy const & t)
{
    tabular_point res;
    for (unsigned int i = 0; i < t.size(); i++)
        res.emplace_back(point {.number = i, .x = t[i].proba, .y = t[i].time});
    return res;
}

tabular_strategy convert_tab_strategy_to_tab_point(tabular_point const & t,
                                                   tabular_strategy const & init)
{
    tabular_strategy res;
    for (auto const & p: t)
        res.push_back(init[p.number]);
    return res;
}

tabular_strategy convex_hull_strategy(tabular_strategy const & t)
{
    return convert_tab_strategy_to_tab_point(
        convex_hull(convert_tab_point_to_tab_strategy(t)), t);
}
