#include "cado.h" // IWYU pragma: keep

#include <cstdlib>

#include "fmt/base.h"

#include "fm.hpp"
#include "generate_factoring_method.hpp"
#include "gmp_aux.h"
#include "tab_fm.hpp"

// coverity[root_function]
int main()
{
    cxx_gmp_randstate state;
    gmp_randseed_ui(state, 42);

    // We exercise the code, without really checking that it works
    tabular_fm tab;
    tab.push_back(factoring_method::from_fields(PM1_METHOD, 0, 20, 100));

    bench_proba(state, tab, 20, 20, 5);
    bench_time(state, tab, 100);

    for (size_t i = 0; i < tab.size(); i++) {
        fmt::print("method {}:", i);
        for (double const t: tab[i].time)
            fmt::print(" {:3.2f}", t);
        fmt::print("\n");
    }

    return EXIT_SUCCESS;
}
