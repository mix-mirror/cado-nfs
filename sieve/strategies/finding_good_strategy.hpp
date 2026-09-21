#ifndef FINDING_GOOD_STRATEGY_HPP
#define FINDING_GOOD_STRATEGY_HPP

#include <istream>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include "generate_strategies.hpp" // strategy_matrix
#include "strategy.hpp"
#include "tab_strategy.hpp"

/* How many times each pair of cofactor sizes was seen, as las reports it
 * with -stats-cofact. */
using cofactor_distribution = std::vector<std::vector<unsigned long>>;

/* One chosen strategy per pair of sizes, empty where the distribution
 * never saw that pair. */
using best_strategies = std::vector<std::vector<std::optional<strategy_t>>>;

strategy_matrix extract_matrix_strat(std::string const & pathname_st,
                                     unsigned int len_abs, unsigned int len_ord);

cofactor_distribution extract_matrix_C(std::istream & is, unsigned int len_abs,
                                       unsigned int len_ord);

best_strategies compute_best_strategy(strategy_matrix const & matrix_strat,
                                      cofactor_distribution const & distrib_C,
                                      unsigned int len_abs,
                                      unsigned int len_ord, double C0);

// to print our final strategies
void strategy_fprint_design(std::ostream & os, strategy_t const & t);

void fprint_final_strategy(std::ostream & os, best_strategies const & res,
                           unsigned int len_abs, unsigned int len_ord);

#endif /* FINDING_GOOD_STRATEGY_HPP */
