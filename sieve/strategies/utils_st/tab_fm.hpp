#ifndef TAB_FM_HPP
#define TAB_FM_HPP

#include <istream>
#include <ostream>
#include <vector>

#include "fmt/base.h"
#include "fmt/ostream.h"

#include "facul_ecm.h"
#include "facul_method.hpp"

#include "fm.hpp"

using tabular_fm = std::vector<factoring_method>;

/* The methods of one family. curve is ignored unless method is
 * EC_METHOD. */
tabular_fm extract_fm_method(tabular_fm const & t, facul_method_code method,
                             ec_parameterization_t curve);

/* Positive if el1 sorts after el2. Beware: this is not a strict weak
 * ordering -- it never reports two methods as equivalent, and it sums
 * signed differences rather than comparing them. It is exposed only so
 * that a test can pin the order that sort_by_proba produces. */
int fm_compare_by_proba(factoring_method const & el1,
                        factoring_method const & el2);

/* Sort by increasing probability, zero methods first. */
void sort_by_proba(tabular_fm & t);

std::istream & operator>>(std::istream & is, tabular_fm &);
std::ostream & operator<<(std::ostream & os, tabular_fm const &);

namespace fmt
{
template <> struct formatter<tabular_fm> : ostream_formatter {
};
} // namespace fmt

#endif /* TAB_FM_HPP */
