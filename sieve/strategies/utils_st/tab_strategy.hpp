#ifndef TAB_STRATEGY_HPP
#define TAB_STRATEGY_HPP

#include <istream>
#include <ostream>
#include <vector>

#include "fmt/base.h"
#include "fmt/ostream.h"

#include "strategy.hpp"

using tabular_strategy = std::vector<strategy_t>;

std::istream & operator>>(std::istream & is, tabular_strategy &);
std::ostream & operator<<(std::ostream & os, tabular_strategy const &);

namespace fmt
{
template <> struct formatter<tabular_strategy> : ostream_formatter {
};
} // namespace fmt

#endif /* TAB_STRATEGY_HPP */
