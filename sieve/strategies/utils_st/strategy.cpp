#include "cado.h" // IWYU pragma: keep

#include <cctype>

#include <istream>
#include <ostream>
#include <utility>

#include "fmt/base.h"

#include "strategy.hpp"
#include "utils_cxx.hpp"

std::ostream & operator<<(std::ostream & os, strategy_t const & s)
{
    for (size_t i = 0; i < s.tab_fm.size(); i++) {
        auto const & p = s.tab_fm[i].params;
        os << p.method << ' ' << p.parameterization << ' ' << p.B1 << ' '
           << p.B2 << ' ' << s.side_of(i) << '\n';
    }
    os << fmt::format("Probability: {:.10f}\n", s.proba);
    os << fmt::format("Time: {:f}\n", s.time);
    return os;
}

std::istream & operator>>(std::istream & is, strategy_t & s)
{
    strategy_t r;

    for (;;) {
        is >> std::ws;
        if (!is.good())
            return is;
        if (!isdigit(is.peek()))
            break;
        unsigned long m[4];
        for (auto & x: m)
            if (!(is >> x))
                return is;
        int side;
        if (!(is >> side))
            return is;
        r.add_fm(factoring_method::from_fields(m[0], m[1], m[2], m[3]), side);
    }

    if (!(is >> std::ws >> expect("Probability:") >> r.proba))
        return is;
    if (!(is >> std::ws >> expect("Time:") >> r.time))
        return is;

    s = std::move(r);
    return is;
}
