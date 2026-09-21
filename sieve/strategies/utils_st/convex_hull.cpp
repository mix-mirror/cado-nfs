#include "cado.h" // IWYU pragma: keep

#include <cmath>

#include <algorithm>

#include "convex_hull.hpp"
#include "point.hpp"
#include "tab_point.hpp"

/* The hull is walked from the leftmost point to the rightmost one. */

static int cmp_double(double a, double b)
{
    double const diff = a - b;
    double const precision = 1e-10;

    if (diff < precision && diff > -precision)
        return 0;
    if (diff < precision)
        return -1;
    return 1;
}

/* p2 is assumed to be to the right of p1, so the adjacent side is
 * positive and we are really comparing slopes. 'scale' brings the two
 * coordinates to a comparable magnitude; see convex_hull(). */
static double compute_angle(point const & p1, point const & p2, double scale)
{
    double const op = p2.y - p1.y;
    double const adj = (p2.x - p1.x) * scale;

    if (adj == 0)
        return op > 0 ? 90 : 360;
    if (op == 0)
        return adj > 0 ? 0 : 180;
    return atan(op / adj);
}

/* The next point of the hull after pt: among those strictly to its
 * right, the one reached with the smallest slope. -1 when pt is the
 * rightmost point. */
static int select_next_point(tabular_point const & t, point const & pt,
                             double scale)
{
    int res = -1;
    double angle_min = 90;

    for (int i = 0; i < (int)t.size(); i++) {
        auto const & elem = t[i];
        if (cmp_double(pt.x, elem.x) == -1) {
            double const angle = compute_angle(pt, elem, scale);
            if (angle < angle_min) {
                angle_min = angle;
                res = i;
            }
        }
    }
    return res;
}

/* return the index of the leftmost point, or -1 if t is empty */
int search_init_point(tabular_point const & t)
{
    if (t.empty())
        return -1;
    int res = 0;
    point pt_min = t[0];
    for (int i = 0; i < (int)t.size(); i++) {
        auto const & elem = t[i];
        if (cmp_double(pt_min.x, elem.x) > 0 ||
            (cmp_double(pt_min.x, elem.x) == 0 &&
             cmp_double(pt_min.y, elem.y) > 0)) {
            res = i;
            pt_min = elem;
        }
    }
    return res;
}

tabular_point convex_hull(tabular_point const & t)
{
    /* the convex hull of the empty set is the empty set. Returning
     * early also avoids the meaningless scaling computation below. */
    if (t.empty())
        return {};

    double minx = INFINITY;
    double maxx = 0;
    double miny = INFINITY;
    double maxy = 0;
    for (auto const & p: t) {
        minx = std::min(minx, p.x);
        miny = std::min(miny, p.y);
        maxx = std::max(maxx, p.x);
        maxy = std::max(maxy, p.y);
    }
    /* x is a probability and y is a time, so they are not commensurate;
     * without this the slope comparison would be dominated by whichever
     * of the two spans the larger range. */
    double const scaley = log(maxy - miny) / log(10);
    double const scalex = log(maxx - minx) / log(10);
    double const scale = pow(10, scaley - scalex);

    tabular_point convex_hull;
    int index = search_init_point(t);
    while (index != -1) {
        convex_hull.push_back(t[index]);
        index = select_next_point(t, t[index], scale);
    }

    return convex_hull;
}
