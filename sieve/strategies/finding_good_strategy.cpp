#include "cado.h" // IWYU pragma: keep

#include <cfloat>
#include <cmath>

#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

#include "fmt/base.h"
#include "fmt/ostream.h"

#include "facul_ecm.h"
#include "facul_method.hpp"
#include "finding_good_strategy.hpp"
#include "fm.hpp"
#include "macros.h"
#include "utils_cxx.hpp"
#include "strategy.hpp"
#include "tab_strategy.hpp"

// #define STATS
/* if defined then print in a file 'result_strat' the different values of
 * s that are tested. */

static double EPSILON_DBL = LDBL_EPSILON;

/************************************************************************/
/*                    EXTRACT MATRIX STRATREGIES                        */
/************************************************************************/

/* Read the strategies of each pair (r0,r1) from the directory
 * 'pathname_st'. */
strategy_matrix extract_matrix_strat(std::string const & pathname_st,
                                     unsigned int len_abs, unsigned int len_ord)
{
    strategy_matrix matrix(len_abs);

    for (unsigned int r0 = 0; r0 < len_abs; r0++) {
        matrix[r0].resize(len_ord);
        for (unsigned int r1 = 0; r1 < len_ord; r1++) {
            auto const name =
                fmt::format("{}/strategies_{}_{}", pathname_st, r0, r1);
            std::ifstream is(name);
            if (!is || !(is >> matrix[r0][r1]))
                throw cado::error("Impossible to read the file '{}'", name);
        }
    }
    return matrix;
}

/************************************************************************/
/*                         EXTRACT THE SET C                            */
/************************************************************************/

/* The distribution of cofactor sizes, as las writes it. Pairs larger
 * than the matrix are dropped. */
cofactor_distribution extract_matrix_C(std::istream & is, unsigned int len_abs,
                                       unsigned int len_ord)
{
    cofactor_distribution matrix_call(len_abs,
                                      std::vector<unsigned long>(len_ord, 0));

    for (;;) {
        is >> std::ws;
        if (is.eof())
            break;
        unsigned int i, j;
        unsigned long c, unused_s;
        if (!(is >> i >> j >> c >> unused_s))
            throw cado::error("Cannot parse the cofactor distribution");
        if (i < len_abs && j < len_ord)
            matrix_call[i][j] = c;
    }

    return matrix_call;
}

/************************************************************************/
/*          SEARCH THE BEST STRATEGIES TO MAXIMIZE Y/T                  */
/************************************************************************/

static unsigned int subroutine_dicho(tabular_strategy const & tab_strat,
                                     double s, unsigned int ind_min,
                                     unsigned int ind_max, double slope_min,
                                     double slope_max)
{
    if (ind_max - ind_min <= 1)
        return ind_min;

    unsigned int const middle = (ind_max + ind_min) / 2;
    double slope_middle;
    if (middle == 0) {
        slope_middle = INFINITY;
    } else {
        auto const & elem1 = tab_strat[middle - 1];
        auto const & elem2 = tab_strat[middle];
        slope_middle = (elem2.proba - elem1.proba) /
                       (elem2.time - elem1.time) * 1000000;
    }

    if (slope_middle < s)
        return subroutine_dicho(tab_strat, s, ind_min, middle, slope_min,
                                slope_middle);
    return subroutine_dicho(tab_strat, s, middle, ind_max, slope_middle,
                            slope_max);
}

static unsigned int
subroutine_compute_slope_yt_dicho(tabular_strategy const & tab_strat, double s)
{
    if (tab_strat.size() == 1)
        return 0;
    return subroutine_dicho(tab_strat, s, 0, tab_strat.size(), INFINITY, 0);
}

static double compute_slope_yt(strategy_matrix const & matrix_strat,
                               cofactor_distribution const & distrib_C,
                               unsigned int len_abs, unsigned int len_ord,
                               double s, double C0)
{
    double Y = 0, T = C0;

    for (unsigned int r1 = 0; r1 < len_abs; r1++) {
        for (unsigned int r2 = 0; r2 < len_ord; r2++) {
            if (distrib_C[r1][r2] > EPSILON_DBL) {
                unsigned int const index =
                    subroutine_compute_slope_yt_dicho(matrix_strat[r1][r2], s);
                Y += double(distrib_C[r1][r2]) *
                     matrix_strat[r1][r2][index].proba;
                T += double(distrib_C[r1][r2]) *
                     matrix_strat[r1][r2][index].time;
            }
        }
    }
    return Y / T;
}

static double sampling_function_interval(
    strategy_matrix const & matrix_strat,
    cofactor_distribution const & distrib_C, unsigned int len_abs,
    unsigned int len_ord, double C0, double init_s, double maxi_s, double pas)
{
#ifdef STATS
    FILE * result_file = fopen("result_strat", "w+");
#endif
    double max_s = 0, max_yt = 0;
    for (double s = init_s; s < maxi_s; s += pas) {
        double const yt =
            compute_slope_yt(matrix_strat, distrib_C, len_abs, len_ord, s, C0);
        if (yt > max_yt) {
            max_s = s;
            max_yt = yt;
        }
#ifdef STATS
        fprintf(result_file, "%lf \t %1.10lf\n", s, yt);
#endif
    }
#ifdef STATS
    fclose(result_file);
#endif
    return max_s;
}

static double sampling_function(strategy_matrix const & matrix_strat,
                                cofactor_distribution const & distrib_C,
                                unsigned int len_abs, unsigned int len_ord,
                                double C0, double init_s, double pas)
{
#ifdef STATS
    FILE * result_file = fopen("result_strat", "w+");
#endif
    double max_s = 0, max_yt = 0;
    double s = init_s;
    int chronos = 0;
    while (chronos < 100) {
        double const yt =
            compute_slope_yt(matrix_strat, distrib_C, len_abs, len_ord, s, C0);
        if (yt > max_yt) {
            chronos = 0;
            max_s = s;
            max_yt = yt;
        } else {
            chronos++;
        }
        s += pas;
#ifdef STATS
        fprintf(result_file, "%lf \t %1.10lf\n", s, yt);
#endif
    }
#ifdef STATS
    fclose(result_file);
#endif

    return max_s;
}

/* Look for the best choice of s: the one that maximizes the number of
 * relations per second of the sieving step (cofactoring plus the sieving
 * that produced the cofactors). */
best_strategies compute_best_strategy(strategy_matrix const & matrix_strat,
                                      cofactor_distribution const & distrib_C,
                                      unsigned int len_abs,
                                      unsigned int len_ord, double C0)
{
    double const step_s = 1;
    // find an interesting interval to search the optimal value of s
    double max_s = sampling_function(matrix_strat, distrib_C, len_abs, len_ord,
                                     C0, 0, step_s);

    // and search it
    double const min = (max_s - step_s > 0) ? max_s - step_s : 0;

    max_s = sampling_function_interval(matrix_strat, distrib_C, len_abs,
                                       len_ord, C0, min, min + 2 * step_s,
                                       0.00010);

    double Y = 0, T = C0;
    best_strategies matrix_res(len_abs);
    for (unsigned int r1 = 0; r1 < len_abs; r1++) {
        matrix_res[r1].resize(len_ord);
        for (unsigned int r2 = 0; r2 < len_ord; r2++) {
            if (distrib_C[r1][r2] < EPSILON_DBL)
                continue;
            unsigned int const index =
                subroutine_compute_slope_yt_dicho(matrix_strat[r1][r2], max_s);
            matrix_res[r1][r2] = matrix_strat[r1][r2][index];

            Y += double(distrib_C[r1][r2]) * matrix_strat[r1][r2][index].proba;
            T += double(distrib_C[r1][r2]) * matrix_strat[r1][r2][index].time;
        }
    }
    fmt::print(" Y = {:f} relations, T = {:f} s., yt = {:1.10f} s/rel.\n", Y,
               T / 1000000, Y ? T / (Y * 1000000) : 0);

    return matrix_res;
}

/************************************************************************/
/*                  DESIGN OUR RESULT                                   */
/************************************************************************/

static char const * method_name(facul_method::parameters const & p)
{
    switch (p.method) {
    case PP1_27_METHOD:
        return "PP1-27";
    case PP1_65_METHOD:
        return "PP1-65";
    case PM1_METHOD:
        return "PM1";
    default: // EC_METHOD
        switch (p.parameterization) {
        case BRENT12:
            return "ECM-B12";
        case MONTY12:
            return "ECM-M12";
        default: // MONTY16
            return "ECM-M16";
        }
    }
}

/* Print the strategies chosen for each pair (r0,r1) that las actually
 * met in the distribution of cofactors. */
void strategy_fprint_design(std::ostream & os, strategy_t const & t)
{
    for (size_t i = 0; i < t.tab_fm.size(); i++) {
        auto const & p = t.tab_fm[i].params;
        fmt::print(os, "[S{}: {}, {}, {} ] ", t.side_of(i), method_name(p),
                   p.B1, p.B2);
    }
    os << "\n";
}

void fprint_final_strategy(std::ostream & os, best_strategies const & res,
                           unsigned int len_abs, unsigned int len_ord)
{
    for (unsigned int r0 = 0; r0 < len_abs; r0++) {
        for (unsigned int r1 = 0; r1 < len_ord; r1++) {
            if (!res[r0][r1])
                continue;
            fmt::print(os, "[r0={}, r1={}] : (p = {:f}, t = {:f})\n", r0, r1,
                       res[r0][r1]->proba, res[r0][r1]->time);
            strategy_fprint_design(os, *res[r0][r1]);
        }
    }
}
