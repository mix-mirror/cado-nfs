#include "cado.h" // IWYU pragma: keep

#include <cstdlib>

#include <sstream>
#include <string>
#include <vector>

#include "fmt/base.h"

#include "fm.hpp"
#include "strategy.hpp"
#include "tab_fm.hpp"
#include "tab_strategy.hpp"

static factoring_method make(unsigned long method, unsigned long curve,
                             unsigned long B1, unsigned long B2)
{
    factoring_method fm = factoring_method::from_fields(method, curve, B1, B2);
    fm.len_p_min = 20;
    fm.proba = {0.5};
    fm.time = {1, 2, 3, 4};
    return fm;
}

/* The file format keeps only the method, the curve, the bounds and the
 * side, so that is all a round trip can preserve. */
static bool same_chain(strategy_t const & a, strategy_t const & b)
{
    if (a.tab_fm.size() != b.tab_fm.size())
        return false;
    for (size_t i = 0; i < a.tab_fm.size(); i++) {
        if (!a.tab_fm[i].same_method_as(b.tab_fm[i]))
            return false;
        if (a.side_of(i) != b.side_of(i))
            return false;
    }
    return true;
}

int main()
{
    int rc = EXIT_SUCCESS;
    auto fail = [&rc](std::string const & what) {
        fmt::print(stderr, "error: {}\n", what);
        rc = EXIT_FAILURE;
    };

    /* a strategy with no side decided yet prints side 0 */
    {
        strategy_t s;
        s.add_fm(make(PM1_METHOD, 0, 315, 2205));
        s.add_fm(make(EC_METHOD, MONTY12, 105, 3255));
        if (s.has_sides())
            fail("a strategy built without sides claims to have them");
        if (s.side_of(0) != 0 || s.side_of(1) != 0)
            fail("a strategy without sides did not default to side 0");
    }

    /* adding with a side keeps the two arrays in step */
    {
        strategy_t s;
        s.add_fm(make(PM1_METHOD, 0, 315, 2205), 1);
        s.add_fm(make(EC_METHOD, MONTY12, 105, 3255), 0);
        if (!s.has_sides())
            fail("a strategy built with sides does not have them");
        if (s.side_of(0) != 1 || s.side_of(1) != 0)
            fail("the sides came back wrong");
    }

    /* mixing the two: the methods added without a side stay on side 0 */
    {
        strategy_t s;
        s.add_fm(make(PM1_METHOD, 0, 315, 2205));
        s.add_fm(make(EC_METHOD, MONTY12, 105, 3255), 1);
        if (!s.has_sides())
            fail("the side array was not filled in to match");
        if (s.side_of(0) != 0 || s.side_of(1) != 1)
            fail("filling in the side array used the wrong value");
    }

    /* round trip through the stream operators */
    {
        tabular_strategy t;

        strategy_t s1;
        s1.add_fm(make(PM1_METHOD, 0, 315, 2205), 0);
        s1.add_fm(make(EC_METHOD, MONTY12, 105, 3255), 1);
        s1.proba = 0.75;
        s1.time = 123.5;
        t.push_back(s1);

        strategy_t s2;
        s2.add_fm(make(EC_METHOD, MONTY16, 0, 0), 1);
        s2.proba = 0;
        s2.time = 0;
        t.push_back(s2);

        std::ostringstream os;
        os << t;

        tabular_strategy back;
        std::istringstream is(os.str());
        if (!(is >> back))
            fail("could not read back what we just wrote");
        if (back.size() != t.size()) {
            fail("the round trip changed the number of strategies");
        } else {
            for (size_t i = 0; i < t.size(); i++) {
                if (!same_chain(t[i], back[i]))
                    fail("the round trip changed a chain of methods");
                /* Probability is written with ten decimals, time with
                 * six; both of these survive exactly. */
                if (back[i].proba != t[i].proba)
                    fail("the round trip changed a probability");
                if (back[i].time != t[i].time)
                    fail("the round trip changed a time");
            }
        }
    }

    /* malformed input fails the stream rather than killing the process */
    {
        std::istringstream is("1 0 315 2205 0\nProbability: oops\n");
        tabular_strategy t;
        if (is >> t)
            fail("a malformed strategy was accepted");
    }

    /* an empty input is an empty table, not a failure */
    {
        std::istringstream is("");
        tabular_strategy t;
        if (!(is >> t) || !t.empty())
            fail("an empty input was not read as an empty table");
    }

    return rc;
}
