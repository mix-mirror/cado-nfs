#ifndef STRATEGY_HPP
#define STRATEGY_HPP

#include <cstddef>

#include <istream>
#include <ostream>
#include <vector>

#include "fmt/base.h"
#include "fmt/ostream.h"

#include "fm.hpp"
#include "tab_fm.hpp"

/* A chain of factoring methods, with the probability that it finds a
 * non-trivial factor and the average time it takes to do so.
 *
 * TODO: this is very nearly a std::vector<facul_method::parameters_with_side>
 * carrying a (probability, time) pair.
 */
struct strategy_t {
    tabular_fm tab_fm;
    double proba = 0;
    double time = 0;

    /* side[i] is the side tab_fm[i] applies to. It stays empty while the
     * strategy is one-sided and no side has been picked yet; printing
     * then falls back to side 0. */
    std::vector<int> side;

    void add_fm(factoring_method const & fm) { tab_fm.push_back(fm); }

    void add_fm(factoring_method const & fm, int s)
    {
        tab_fm.push_back(fm);
        /* Methods may have been added without a side before this one;
         * they are taken to be on side 0, as they always have been. */
        side.resize(tab_fm.size(), 0);
        side.back() = s;
    }

    bool has_sides() const { return side.size() == tab_fm.size(); }
    int side_of(size_t i) const { return has_sides() ? side[i] : 0; }
};

std::istream & operator>>(std::istream & is, strategy_t &);
std::ostream & operator<<(std::ostream & os, strategy_t const &);

namespace fmt
{
template <> struct formatter<strategy_t> : ostream_formatter {
};
} // namespace fmt

#endif /* STRATEGY_HPP */
