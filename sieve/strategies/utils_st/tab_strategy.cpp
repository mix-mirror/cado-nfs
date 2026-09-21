#include "cado.h" // IWYU pragma: keep

#include <istream>
#include <ostream>
#include <utility>

#include "tab_strategy.hpp"

std::ostream & operator<<(std::ostream & os, tabular_strategy const & t)
{
    for (auto const & s: t)
        os << s;
    return os;
}

std::istream & operator>>(std::istream & is, tabular_strategy & t)
{
    for (;;) {
        is >> std::ws;
        if (is.eof())
            break;
        if (!is.good())
            return is;
        strategy_t s;
        if (!(is >> s))
            return is;
        t.push_back(std::move(s));
    }
    return is;
}
