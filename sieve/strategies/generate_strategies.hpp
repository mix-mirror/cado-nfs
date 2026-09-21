#ifndef GENERATE_STRATEGIES_HPP
#define GENERATE_STRATEGIES_HPP

#include <string>
#include <vector>

#include "decomp.hpp"
#include "fm.hpp"
#include "strategy.hpp"
#include "tab_decomp.hpp"
#include "tab_fm.hpp"
#include "tab_point.hpp"
#include "tab_strategy.hpp"

/************************************************************************/
/*                      COLLECT DATA FOR ONLY ONE COFACTOR              */
/************************************************************************/

double compute_proba_method_one_decomp(decomp const & dec,
                                       factoring_method const & fm);

double compute_proba_strategy(tabular_decomp const & init_tab,
                              strategy_t const & strat,
                              unsigned int len_p_min, unsigned int len_p_max);

double compute_time_strategy(tabular_decomp const & init_tab,
                             strategy_t const & strat, unsigned int r);

/************************************************************************/
/*                   GENERATE MATRIX                                    */
/************************************************************************/

tabular_strategy generate_strategies_oneside(tabular_decomp const & tab_decomp,
                                             factoring_method const & zero,
                                             tabular_fm const & pm1,
                                             tabular_fm const & pp1,
                                             tabular_fm const & ecm,
                                             int nb_curve, unsigned long lim,
                                             unsigned int lpb, unsigned int r);

tabular_strategy generate_strategy_r0_r1(tabular_strategy const & strat_r0,
                                         tabular_strategy const & strat_r1);

/* matrix[r0][r1], with 0 <= r0 <= mfb0 and 0 <= r1 <= mfb1 */
using strategy_matrix = std::vector<std::vector<tabular_strategy>>;

strategy_matrix generate_matrix(std::string const & name_directory_decomp,
                                tabular_fm const & pm1, tabular_fm const & pp1,
                                tabular_fm const & ecm, int nb_curve,
                                unsigned long lim0, unsigned int lpb0,
                                unsigned int mfb0, unsigned long lim1,
                                unsigned int lpb1, unsigned int mfb1);

/************************************************************************/
/*                      CONVEX_HULL_ST                                  */
/************************************************************************/

tabular_point convert_tab_point_to_tab_strategy(tabular_strategy const & t);

tabular_strategy convert_tab_strategy_to_tab_point(tabular_point const & t,
                                                   tabular_strategy const & init);

tabular_strategy convex_hull_strategy(tabular_strategy const & t);

#endif /* GENERATE_STRATEGIES_HPP */
