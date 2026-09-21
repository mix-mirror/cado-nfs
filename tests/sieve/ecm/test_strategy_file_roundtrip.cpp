#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <sstream>
#include <string>
#include <exception>
#include <vector>

#include "fmt/base.h"
#include "fmt/format.h"

#include "ecm/facul_ecm.h"
#include "ecm/facul_method.hpp"
#include "ecm/facul_strategies.hpp"
#include "macros.h"

/* The strategy file is written by finalst, in sieve/strategies, and read
 * back here. Until now only the reading half was tested. This writes one
 * method of every family the parser knows about, reads the result back,
 * and checks that nothing changed on the way.
 */

static std::vector<facul_method::parameters_with_side> every_family()
{
    std::vector<facul_method::parameters_with_side> v;
    v.emplace_back(0, PM1_METHOD, 315UL, 2205UL);
    v.emplace_back(1, PP1_27_METHOD, 525UL, 3255UL);
    v.emplace_back(0, PP1_65_METHOD, 975UL, 107250UL);
    for (auto const para: {BRENT12, MONTY12, MONTY16, MONTYTWED12})
        v.emplace_back(1, EC_METHOD, 105UL, 3255UL, para, 2UL, 1);
    return v;
}

int main()
{
    int rc = EXIT_SUCCESS;
    auto fail = [&rc](std::string const & what) {
        fmt::print(stderr, "error: {}\n", what);
        rc = EXIT_FAILURE;
    };

    auto const methods = every_family();

    /* every family must have a name, and they must all differ */
    for (size_t i = 0; i < methods.size(); i++)
        for (size_t j = i + 1; j < methods.size(); j++)
            if (facul_method_name(methods[i]) == facul_method_name(methods[j]))
                fail(fmt::format("two families share the name {}",
                                 facul_method_name(methods[i])));

    /* write a one-entry strategy file the way finalst does */
    std::ostringstream os;
    os << "[r0=30, r1=30] : (p = 0.500000, t = 1.000000)\n";
    for (auto const & m: methods)
        os << m;
    os << "\n";

    /* and read it back through the real parser */
    std::string const text = os.str();

    // coverity[secure_temp]
    FILE * f = tmpfile();
    DIE_ERRNO_DIAG(f == nullptr, "tmpfile(%s)", "");
    fputs(text.c_str(), f);
    fflush(f);

    std::vector<unsigned long> const B(2, 1UL << 20);
    std::vector<unsigned int> const lpb(2, 25);
    std::vector<unsigned int> const mfb(2, 50);

    try {
        auto const F = facul_strategies(B, lpb, mfb, true, f, 0);
        std::vector<unsigned int> const index {30, 30};
        auto const & chain = F(index);
        if (chain.size() != methods.size()) {
            fail(fmt::format("wrote {} methods, read back {}", methods.size(),
                             chain.size()));
        } else {
            for (size_t i = 0; i < methods.size(); i++) {
                auto const & w = methods[i];
                auto const & r = chain[i];
                if (r.side != w.side)
                    fail(fmt::format("method {}: side {} became {}", i, w.side,
                                     r.side));
                if (r.method->method != w.method)
                    fail(fmt::format("method {}: {} did not survive the round trip",
                                     i, facul_method_name(w)));
            }
        }
    } catch (std::exception const & e) {
        fail(fmt::format("the parser rejected what we wrote: {}\n{}", e.what(),
                         text));
    }

    fclose(f);
    return rc;
}
