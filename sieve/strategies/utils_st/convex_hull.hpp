#ifndef CONVEX_HULL_HPP
#define CONVEX_HULL_HPP

#include <vector>

#include "point.hpp"
#include "tab_point.hpp"

/* The lower convex hull of a set of (x, y) points, walked from the
 * leftmost point rightwards.
 *
 * Everything in this directory that has to choose between candidates
 * trading a probability x off against a cost y goes through this: it is
 * the frontier of the best trade-offs, and nothing below it exists. */
tabular_point convex_hull(tabular_point const & t);

/* the leftmost point of t, or -1 if t is empty */
int search_init_point(tabular_point const & t);

/* Reduce a collection to the elements whose (x, y) lie on that hull,
 * keeping them in hull order.
 *
 * This is the shape every caller wants: they hold strategies, or
 * factoring methods, and want the ones on the trade-off frontier. The
 * point type exists only to carry an index back to the real object, so
 * it stays inside here. */
template <typename T, typename XY>
std::vector<T> pareto_reduce(std::vector<T> const & v, XY xy)
{
    tabular_point points;
    points.reserve(v.size());
    for (unsigned int i = 0; i < v.size(); i++) {
        auto const [x, y] = xy(v[i]);
        points.emplace_back(point {.number = i, .x = x, .y = y});
    }

    std::vector<T> res;
    for (auto const & p: convex_hull(points))
        res.push_back(v[p.number]);
    return res;
}

#endif /* CONVEX_HULL_HPP */
