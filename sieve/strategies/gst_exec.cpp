#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdlib>

#include <exception>
#include <fstream>
#include <string>
#include <vector>

#include "fmt/base.h"

#include "facul_ecm.h"
#include "facul_method.hpp"
#include "fm.hpp"
#include "gen_decomp.hpp"
#include "generate_strategies.hpp"
#include "macros.h"
#include "params.hpp"
#include "tab_decomp.hpp"
#include "cado_main.hpp"
#include "tab_fm.hpp"
#include "tab_strategy.hpp"

/************************************************************************/
/*                            USAGE                                     */
/************************************************************************/

static void declare_usage(cxx_param_list & pl)
{
    pl.declare_usage_header("This binary allows to build the best strategies for each couple "
            "(r0, r1)\n"
            "where (r0,r1) are the bits size for our couple of cofactors.\n");

    pl.declare_usage("gdc",
                          "to precompute all decompositions of cofactors of "
                          "mfb bits given that \n "
                          "\t \t it has no prime divisors less than lim. So, "
                          "you must specify these options:\n"
                          "\t \t -lim0, -mfb0\n");
    pl.declare_usage("gst_r",
                          "to precompute the best strategies \n"
                          "\t \t for one bit size cofactor.\n "
                          "\t \t You must specify these options:\n"
                          "\t \t -lim0, -lpb0, -r0, -ncurves, -decomp\n");
    pl.declare_usage("gst",
        "to merge two (or all) precomputing did by the option 'gst_r',\n "
        "\t \t and thus find the best strategie(s) for one (or each) couple "
        "(r0,r1).\n"
        "\t \t So, you must specify these options:\n"
        "\t \t -lim0, -lim1 ,-in, and ((-r0 -r1) or (-mfb0, -mfb1))\n");

    pl.declare_usage("ncurves", "controls number of curves.\n");
    pl.declare_usage("r0",
                          "set the bit size of the studied cofactor to r0.\n");
    pl.declare_usage("r1",
                          "set the bit size of the studied cofactor to r1.\n");
    pl.declare_usage("lim0",
                          "set rationnal factor base bound to lim0\n");
    pl.declare_usage("lim1",
                          "set algebraic factor base bound to lim1\n");
    pl.declare_usage("lpb0",
                          "set rational large prime bound to 2^lpb0");
    pl.declare_usage("lpb1",
                          "set algebraic large prime bound to 2^lpb1");
    pl.declare_usage("mfb0", "set the first cofactor bound to 2^mfb0");
    pl.declare_usage("mfb1",
                          "set the second cofactor bound to 2^mfb1");
    pl.declare_usage("in",
        "to locate the file which contains\n "
        "\t \t our factoring methods, or locate the directory \n"
        "\t \t where the precomputed files for option 'gst' are stored");
    pl.declare_usage("out",
                          "to locate the directory where the "
                          "file(s) will be stored.");
    pl.declare_usage("decomp",
        "to locate the file or the directory , according to\n"
        "\t \t if you need one or several files,\n"
        "\t \t which contain(s) the file(s) of cofactors decompositions.");
}

/************************************************************************/
/*     MAIN                                                             */
/************************************************************************/

namespace
{

/* Everything gst reads and writes goes through these, so that the file
 * names stay in one place. */
std::string oneside_filename(std::string const & dir, unsigned long lim, int r)
{
    return fmt::format("{}/strategies{}_{}", dir, lim, r);
}

std::string pair_filename(std::string const & dir, int r0, int r1)
{
    return fmt::format("{}/strategies_{}_{}", dir, r0, r1);
}

template <typename T> T read_from(std::string const & filename)
{
    std::ifstream is(filename);
    T res;
    if (!is || !(is >> res))
        throw cado::error("Cannot read {}", filename);
    return res;
}

template <typename T> void write_to(std::string const & filename, T const & x)
{
    std::ofstream os(filename);
    if (!os || !(os << x))
        throw cado::error("Cannot write {}", filename);
}

/* The six families gst knows about, extracted from one methods file.
 * They are what the strategy generator chains together. */
struct method_families {
    tabular_fm pm1;
    tabular_fm pp1;
    tabular_fm ecm;

    method_families(cxx_param_list & pl, std::string const & filename)
    {
        auto const all = read_from<tabular_fm>(filename);

        auto const pp1_27 = extract_fm_method(all, PP1_27_METHOD, BRENT12);
        auto const pp1_65 = extract_fm_method(all, PP1_65_METHOD, BRENT12);
        auto const ecm_m16 = extract_fm_method(all, EC_METHOD, MONTY16);
        auto const ecm_m12 = extract_fm_method(all, EC_METHOD, MONTY12);
        auto const ecm_b12 = extract_fm_method(all, EC_METHOD, BRENT12);

        pm1 = extract_fm_method(all, PM1_METHOD, BRENT12);
        pp1 = pp1_27;
        pp1.insert(pp1.end(), pp1_65.begin(), pp1_65.end());

        /* Each family is indexed at [0] below and then fed to the convex
         * hull. A family the input file does not mention at all is a user
         * error, not something we can paper over. */
        struct {
            tabular_fm const & data;
            char const * name;
        } const check[] = {
            {pm1, "PM1"},
            {pp1, "PP1-27 or PP1-65"},
            {ecm_b12, "ECM-B12"},
            {ecm_m12, "ECM-M12"},
            {ecm_m16, "ECM-M16"},
        };
        for (auto const & f: check)
            if (f.data.empty())
                pl.fail("Error: file {} contains no {} method\n", filename,
                        f.name);

        /* Every family needs a zero method, so that the generator can
         * leave a slot empty. */
        auto const zero = factoring_method::from_fields(0, 0, 0, 0);
        auto with_zero = [&zero](tabular_fm t) {
            if (t[0].params.B1 != 0)
                t.push_back(zero);
            return t;
        };

        pm1 = with_zero(std::move(pm1));
        pp1 = with_zero(std::move(pp1));

        ecm = with_zero(ecm_m16);
        auto const b12 = with_zero(ecm_b12);
        auto const m12 = with_zero(ecm_m12);
        ecm.insert(ecm.end(), b12.begin(), b12.end());
        ecm.insert(ecm.end(), m12.begin(), m12.end());

        sort_by_proba(pm1);
        sort_by_proba(pp1);
        sort_by_proba(ecm);
    }
};

/* -gdc: all the decompositions of an mfb-bit cofactor with no prime
 * factor below lim. */
void subcommand_gdc(cxx_param_list & pl, unsigned long lim0, int mfb0)
{
    if (lim0 == 0)
        pl.fail("Error: parameter -lim0 is mandatory\n");
    if (mfb0 == -1)
        pl.fail("Error: parameter -mfb0 is mandatory\n");

    auto const res = generate_all_decomp(mfb0, lim0);

    char const * file_out = pl.lookup_old("out");
    if (file_out)
        write_to(file_out, res);
    else
        fmt::print("{}", res);
}

/* -gst_r: the best strategies for one cofactor size on one side. */
void subcommand_gst_r(cxx_param_list & pl, std::string const & dir_out,
                      method_families const & fam, unsigned long lim0,
                      int lpb0, int ncurves)
{
    int r0 = -1;
    pl.parse("r0", r0);
    if (r0 == -1 || lpb0 == -1 || lim0 == 0 || ncurves == -1)
        pl.fail("Error: parameters -r0 -lim0 -lpb0 -ncurves are mandatories.\n");

    unsigned int const fbb0 = ceil(log2((double)(lim0 + 1)));
    tabular_decomp tab_decomp;
    if ((unsigned int)r0 >= 2 * fbb0 - 1)
        tab_decomp = read_from<tabular_decomp>(
            pl.parse_mandatory<std::string>("decomp"));

    auto const zero = factoring_method::from_fields(0, 0, 0, 0);
    auto const res =
        generate_strategies_oneside(tab_decomp, zero, fam.pm1, fam.pp1,
                                    fam.ecm, ncurves, lim0, lpb0, r0);

    write_to(oneside_filename(dir_out, lim0, r0), res);
}

/* -gst: merge the one-sided precomputations of the two sides. */
void subcommand_gst(cxx_param_list & pl, std::string const & dir_out,
                    unsigned long lim0, unsigned long lim1, int mfb0, int mfb1)
{
    char const * dir_in = pl.lookup_old("in");
    if (dir_in == nullptr)
        pl.fail("Error: parameter -in is mandatory\n");
    if (lim0 == 0)
        pl.fail("Error: parameter -lim0 is mandatory\n");
    if (lim1 == 0)
        pl.fail("Error: parameter -lim1 is mandatory\n");

    int r0 = -1, r1 = -1;
    pl.parse("r0", r0);
    pl.parse("r1", r1);

    if (r0 != -1 && r1 != -1) {
        // just one pair
        auto const strat_r0 =
            read_from<tabular_strategy>(oneside_filename(dir_in, lim0, r0));
        auto const strat_r1 =
            read_from<tabular_strategy>(oneside_filename(dir_in, lim1, r1));
        write_to(pair_filename(dir_out, r0, r1),
                 generate_strategy_r0_r1(strat_r0, strat_r1));
        return;
    }

    if (mfb0 == -1)
        pl.fail("Error: parameter -mfb0 is mandatory\n");
    if (mfb1 == -1)
        pl.fail("Error: parameter -mfb1 is mandatory\n");

    std::vector<tabular_strategy> data_rat(mfb0 + 1);
    for (int r = 0; r <= mfb0; r++)
        data_rat[r] =
            read_from<tabular_strategy>(oneside_filename(dir_in, lim0, r));

    for (int b = 0; b <= mfb1; b++) {
        auto const strat_r1 =
            read_from<tabular_strategy>(oneside_filename(dir_in, lim1, b));
        for (int a = 0; a <= mfb0; a++)
            write_to(pair_filename(dir_out, a, b),
                     generate_strategy_r0_r1(data_rat[a], strat_r1));
    }
}

/* No switch: do the whole matrix without any precomputation. */
void subcommand_matrix(cxx_param_list & pl, std::string const & dir_out,
                       method_families const & fam, unsigned long lim0,
                       int lpb0, int mfb0, unsigned long lim1, int lpb1,
                       int mfb1, int ncurves)
{
    if (lim0 == 0 || lpb0 == -1 || mfb0 == -1 || lim1 == 0 || lpb1 == -1 ||
        mfb1 == -1 || ncurves == -1)
        pl.fail("Error: parameters -lim0 -lpb0 -mfb0"
                " -lim1 -lpb1 -mfb1 -ncurves are mandatories\n");

    char const * decomp = pl.lookup_old("decomp");
    if (decomp == nullptr)
        pl.fail("Error: parameter -decomp is mandatory\n");

    auto const matrix =
        generate_matrix(decomp, fam.pm1, fam.pp1, fam.ecm, ncurves, lim0, lpb0,
                        mfb0, lim1, lpb1, mfb1);

    for (int r0 = 0; r0 <= mfb0; r0++)
        for (int r1 = 0; r1 <= mfb1; r1++)
            write_to(pair_filename(dir_out, r0, r1), matrix[r0][r1]);
}

} // namespace

static int gst(int argc, char const * argv[])
{
    cxx_param_list pl;
    declare_usage(pl);

    pl.configure_switch("gdc");
    pl.configure_switch("gst");
    pl.configure_switch("gst_r");

    pl.process_command_line_and_extra_parameter_files(argc, argv);

    unsigned long lim0 = 0;
    unsigned long lim1 = 0;
    int mfb0 = -1, mfb1 = -1;
    int lpb0 = -1, lpb1 = -1;
    int ncurves = -1;

    pl.parse("lim0", lim0);
    pl.parse("lim1", lim1);
    pl.parse("lpb0", lpb0);
    pl.parse("lpb1", lpb1);
    pl.parse("mfb0", mfb0);
    pl.parse("mfb1", mfb1);
    pl.parse("ncurves", ncurves);

    if (pl.parse<int>("-gdc")) {
        subcommand_gdc(pl, lim0, mfb0);
        return EXIT_SUCCESS;
    }

    char const * out = pl.lookup_old("out");
    std::string const dir_out = out ? out : "./";

    if (pl.parse<int>("-gst")) {
        subcommand_gst(pl, dir_out, lim0, lim1, mfb0, mfb1);
        return EXIT_SUCCESS;
    }

    char const * name_file_in = pl.lookup_old("in");
    if (name_file_in == nullptr)
        pl.fail("Error: parameter -in is mandatory\n");
    method_families const fam(pl, name_file_in);

    if (pl.parse<int>("-gst_r"))
        subcommand_gst_r(pl, dir_out, fam, lim0, lpb0, ncurves);
    else
        subcommand_matrix(pl, dir_out, fam, lim0, lpb0, mfb0, lim1, lpb1, mfb1,
                          ncurves);

    return EXIT_SUCCESS;
}

/* Bad parameters and bad input files are reported by an exception;
 * cado::main_wrapper is what turns it into a message on stderr and a
 * nonzero exit status, rather than whatever the C++ runtime prints
 * before it aborts -- which differs between libstdc++ and libc++, and
 * can be nothing at all.
 */
// coverity[root_function]
int main(int argc, char const * argv[])
{
    return cado::main_wrapper(gst, argc, argv);
}
