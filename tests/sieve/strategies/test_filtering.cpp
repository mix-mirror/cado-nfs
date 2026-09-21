#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdlib>

#include "fmt/base.h"

#include "fm.hpp"
#include "generate_factoring_method.hpp"
#include "gmp_aux.h"
#include "tab_fm.hpp"

/* check that the spread of the remaining methods is homogeneous */
static bool check_filt(tabular_fm const & res, unsigned int init_nb_method)
{
    size_t const len = res.size();
    if (len <= 2) // fewer than two points is not representative
        return true;
    double aver_gap_prat = 0;
    double const aver_gap_theo = init_nb_method / double(len);
    for (size_t i = 0; i + 1 < len; i++)
        aver_gap_prat +=
            double(res[i + 1].params.method) - double(res[i].params.method);

    aver_gap_prat /= double(len - 1);
    double perc_gap = (aver_gap_prat - aver_gap_theo) / aver_gap_theo;
    if (perc_gap < 0)
        perc_gap *= -1;
    return perc_gap <= 0.4;
}

// coverity[root_function]
int main()
{
    cxx_gmp_randstate state;

    unsigned int const nb_fm = 10;
    unsigned int const final_nb_fm = 4;

    tabular_fm tab;
    for (unsigned int i = 0; i < nb_fm; i++) {
        factoring_method fm = factoring_method::from_fields(
            i, 0, i * (1 + gmp_urandomm_ui(state, 10)),
            i * (1 + gmp_urandomm_ui(state, 10)));
        fm.len_p_min = 0;
        fm.proba.resize(4);
        for (int j = 0; j < 4; j++)
            fm.proba[j] = i / (double(nb_fm) + 1) + 0.01 * j * i;
        fm.time.resize(4);
        for (int j = 0; j < 4; j++)
            fm.time[j] = double(i * i * i) * pow(10, j);
        tab.push_back(std::move(fm));
    }

    tabular_fm const res = filtering(tab, final_nb_fm);

    if (!check_filt(res, nb_fm)) {
        fmt::print(stderr, "the remaining methods are not spread out\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
