#include "cado.h" // IWYU pragma: keep

#include <cstdlib>

#include <sstream>
#include <string>

#include "fmt/base.h"

#include "fm.hpp"
#include "gmp_aux.h"
#include "macros.h"
#include "random_distributions.hpp"
#include "tab_fm.hpp"

static factoring_method make(unsigned long method, unsigned long curve,
                             unsigned long B1, unsigned long B2,
                             std::vector<double> proba,
                             std::vector<double> time)
{
    factoring_method fm =
        factoring_method::from_fields(method, curve, B1, B2);
    fm.len_p_min = 20;
    fm.proba = std::move(proba);
    fm.time = std::move(time);
    return fm;
}

static bool same(factoring_method const & a, factoring_method const & b)
{
    return a.same_method_as(b) && a.len_p_min == b.len_p_min &&
           a.proba == b.proba && a.time == b.time;
}

/* is the table sorted the way sort_by_proba() promises? */
static bool is_sorted_by_proba(tabular_fm const & t)
{
    for (size_t i = 0; i + 1 < t.size(); i++)
        if (fm_compare_by_proba(t[i], t[i + 1]) > 0)
            return false;
    return true;
}

int main()
{
    cxx_gmp_randstate state;
    int rc = EXIT_SUCCESS;
    auto fail = [&rc](std::string const & what) {
        fmt::print(stderr, "error: {}\n", what);
        rc = EXIT_FAILURE;
    };

    /* the four numbers of the file format land where they should */
    {
        auto const pm1 = make(PM1_METHOD, 0, 315, 2205, {0.5}, {1, 2, 3, 4});
        if (pm1.method() != PM1_METHOD || pm1.B1() != 315 || pm1.B2() != 2205)
            fail("the method fields did not survive from_fields()");
        auto const ecm = make(EC_METHOD, MONTY12, 105, 3255, {0.5}, {1});
        if (ecm.curve() != MONTY12)
            fail("the curve did not survive from_fields()");
        /* the curve slot is meaningless for PM1 but is carried anyway */
        auto const odd = make(PM1_METHOD, MONTY12, 315, 2205, {0.5}, {1});
        if (odd.params.parameterization != MONTY12)
            fail("the curve slot of a non-ECM method was dropped");
    }

    /* the zero method */
    {
        auto z = make(EC_METHOD, MONTY12, 105, 3255, {0.5, 0.25}, {1, 2});
        if (z.is_zero())
            fail("a method with B1 != 0 reported itself as zero");
        z.put_zero();
        if (!z.is_zero())
            fail("put_zero() did not produce a zero method");
        for (double p: z.proba)
            if (p != 0)
                fail("put_zero() left a non-zero probability");
    }

    /* time_for() picks the right bucket and clamps to what we have */
    {
        auto const fm = make(PM1_METHOD, 0, 315, 2205, {0.5}, {10, 20, 30, 40});
        if (fm.time_for(1) != 10)
            fail("time_for() missed the first bucket");
        if (fm.time_for(1000) != 40)
            fail("time_for() did not clamp to the last bucket");
    }

    /* round trip through the stream operators */
    {
        tabular_fm t;
        t.push_back(make(PM1_METHOD, 0, 315, 2205, {0.5, 0.25}, {1, 2, 3, 4}));
        t.push_back(make(EC_METHOD, MONTY12, 105, 3255, {0.9}, {5}));
        t.push_back(make(EC_METHOD, MONTY16, 0, 0, {0}, {0}));

        std::ostringstream os;
        os << t;

        tabular_fm back;
        std::istringstream is(os.str());
        if (!(is >> back))
            fail("could not read back what we just wrote");
        if (back.size() != t.size()) {
            fail("the round trip changed the number of methods");
        } else {
            for (size_t i = 0; i < t.size(); i++)
                if (!same(t[i], back[i]))
                    fail("the round trip changed a method");
        }
    }

    /* comments are skipped -- the old fscanf reader could not do this */
    {
        std::istringstream is("# a comment, with 123 digits in it\n"
                              "1 0 315 2205 | 20 0.500000 | 1.000000 |\n"
                              "# another one\n");
        tabular_fm t;
        if (!(is >> t))
            fail("a comment line made the reader fail");
        if (t.size() != 1)
            fail("a comment line was read as a method");
    }

    /* malformed input fails the stream rather than killing the process */
    {
        std::istringstream is("1 0 315 2205 | 20 0.5 ; 1.0 |\n");
        tabular_fm t;
        if (is >> t)
            fail("a malformed record was accepted");
    }

    /* extract_fm_method picks one family */
    {
        tabular_fm t;
        t.push_back(make(PM1_METHOD, 0, 315, 2205, {0.5}, {1}));
        t.push_back(make(EC_METHOD, MONTY12, 105, 3255, {0.9}, {1}));
        t.push_back(make(EC_METHOD, MONTY16, 105, 3255, {0.9}, {1}));
        t.push_back(make(EC_METHOD, MONTY12, 205, 3255, {0.9}, {1}));

        if (extract_fm_method(t, PM1_METHOD, BRENT12).size() != 1)
            fail("extract_fm_method did not find the single PM1");
        if (extract_fm_method(t, EC_METHOD, MONTY12).size() != 2)
            fail("extract_fm_method did not find both MONTY12 curves");
        if (!extract_fm_method(t, PP1_27_METHOD, BRENT12).empty())
            fail("extract_fm_method invented a PP1 method");
    }

    /* sorting */
    {
        tabular_fm t;
        t.push_back(make(PM1_METHOD, 0, 315, 2205, {0.9, 0.8, 0.7, 0.6, 0.3}, {1}));
        t.push_back(make(PM1_METHOD, 0, 315, 2205, {0.8999, 0.8001, 0.6999, 0.55, 0.2999}, {1}));
        t.push_back(make(PM1_METHOD, 0, 315, 2205, {0.5, 0.399, 0.19, 0.0006, 0.00003}, {1}));
        t.push_back(make(PM1_METHOD, 0, 315, 2205, {0.2, 0.18, 0.07, 0.0006, 0.0000003}, {1}));
        t.push_back(make(PM1_METHOD, 0, 315, 2205, {0.6, 0.5, 0.4, 0.3, 0.2}, {1}));
        /* and one random table, to exercise the comparator more widely */
        for (int i = 0; i < 10; i++) {
            std::vector<double> p(5);
            for (auto & x: p)
                x = random_uniform(state);
            t.push_back(make(PM1_METHOD, 0, 315, 2205, std::move(p), {1}));
        }

        sort_by_proba(t);
        if (!is_sorted_by_proba(t))
            fail("sort_by_proba() left the table unsorted");

        if (t.size() > 3) {
            std::swap(t[1], t[3]);
            if (is_sorted_by_proba(t))
                fail("swapping two entries left the table sorted");
        }
    }

    return rc;
}
