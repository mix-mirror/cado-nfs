#include "cado.h" // IWYU pragma: keep

#include <cfloat>
#include <istream>
#include <limits>
#include <ostream>
#include <utility>

#include <algorithm>
#include "tab_fm.hpp"

static double const EPSILON_DBL = 0.000001;

tabular_fm extract_fm_method(tabular_fm const & t, facul_method_code method,
                             ec_parameterization_t curve)
{
    tabular_fm res;
    for (auto const & fm: t) {
        if (fm.params.method != method)
            continue;
        if (method == EC_METHOD && fm.params.parameterization != curve)
            continue;
        res.push_back(fm);
    }
    return res;
}

/************************************************************************/
/*                      SORT                                            */
/************************************************************************/

/* Positive if el1 is "greater" than el2.
 *
 * Note that this is not a strict weak ordering: it never reports two
 * methods as equivalent, and it sums signed differences rather than
 * comparing them. Replacing it, or handing it to std::sort, changes
 * which order comes out. It is kept as it was, along with the selection
 * sort below, so that the files this produces do not move. */
int fm_compare_by_proba(factoring_method const & el1,
                        factoring_method const & el2)
{
    if (el1.is_zero())
        return -1;
    if (el2.is_zero())
        return 1;
    /* assumes that len_p_min is the same for both */
    size_t const len = std::min(el1.proba.size(), el2.proba.size());
    double diff_proba = 0;
    for (size_t i = 0; i < len; i++)
        if (el1.proba[i] > EPSILON_DBL && el2.proba[i] > EPSILON_DBL)
            diff_proba += el1.proba[i] - el2.proba[i];

    return (diff_proba > EPSILON_DBL) ? 1 : -1;
}

void sort_by_proba(tabular_fm & t)
{
    for (size_t max = t.size(); max > 0; max--) {
        size_t index_max = 0;
        for (size_t i = 0; i < max; i++)
            if (fm_compare_by_proba(t[i], t[index_max]) > 0)
                index_max = i;
        std::swap(t[max - 1], t[index_max]);
    }
}

/************************************************************************/
/*                      I/O                                             */
/************************************************************************/

std::ostream & operator<<(std::ostream & os, tabular_fm const & t)
{
    for (auto const & fm: t)
        os << fm;
    return os;
}

std::istream & operator>>(std::istream & is, tabular_fm & t)
{
    for (;;) {
        is >> std::ws;
        if (is.eof())
            break;
        if (!is.good())
            return is;
        if (is.peek() == '#') {
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }
        factoring_method fm;
        if (!(is >> fm))
            return is;
        t.push_back(std::move(fm));
    }
    return is;
}
