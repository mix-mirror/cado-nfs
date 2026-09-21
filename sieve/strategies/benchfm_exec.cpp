#include "cado.h" // IWYU pragma: keep

#include <cstdlib>

#include <fstream>
#include <gmp.h>

#include "tab_fm.hpp"
#include "generate_factoring_method.hpp"
#include "params.hpp"
#include "cado_main.hpp"
#include "utils_cxx.hpp"

static void declare_usage(cxx_param_list & pl)
{
    pl.declare_usage_header("This binary allows to do a bench on the probabilities and/or the times\n"
			    "for each size of prime number from 'lb'.\n");

    pl.declare_usage("lb",
			  "to begin the benchmark with prime numbers of 'lb' bits.");
    pl.declare_usage("p", "to bench the probabilities.");
    pl.declare_usage("t", "to bench the times.");
    pl.declare_usage("N", "number of tests to run for each bench");
    pl.declare_usage("in",
			  "to locate the file which contains our factoring methods.");
    pl.declare_usage("f",
			  "to keep only a number of factoring methods.");
    pl.declare_usage("out",
			  "to locate the file which contains our benchmark.");
    pl.declare_usage("seed", "random seed");

}

/************************************************************************/
/*                      MAIN */
/************************************************************************/

static int main_(int argc, char const * argv[]);

// coverity[root_function]
int main(int argc, char const * argv[])
{
    return cado::main_wrapper(main_, argc, argv);
}

static int main_(int argc, char const * argv[])
{
    int nb_test = 0;
    cxx_param_list pl;
    declare_usage(pl);
    pl.configure_switch("p");
    pl.configure_switch("t");

    pl.process_command_line_and_extra_parameter_files(argc, argv);

    //default values
    int len_p_min = -1; //default_value
    int final_nb_fm = -1; //default value
    
    int const opt_proba = pl.parse<int>("-p");
    int const opt_time = pl.parse<int>("-t");
    pl.parse("lb", len_p_min);
    pl.parse("f", final_nb_fm);
    pl.parse("N", nb_test);

    const char *pathname_in;
    const char *pathname_out;
    if ((pathname_in = pl.lookup_old("in")) == NULL)
        pl.fail("missing argument -in");
    if ((pathname_out = pl.lookup_old("out")) == NULL)
        pl.fail("missing argument -out");

    cxx_gmp_randstate state;
    unsigned long seed = 0;
    if (pl.parse("seed", seed))
        gmp_randseed_ui(state, seed);

    tabular_fm c;
    {
        std::ifstream file_in(pathname_in);
        if (!file_in || !(file_in >> c))
            throw cado::error("impossible to read {}", pathname_in);
    }

    if (opt_proba)
	{
	    if (len_p_min == -1)
                pl.fail("missing argument -lb");
	    bench_proba(state, c, len_p_min, 0, nb_test);
	}
    if (opt_time)
	bench_time(state, c, nb_test);

    if (final_nb_fm != -1)
	c = filtering (c, final_nb_fm);

    std::ofstream file_out(pathname_out);
    if (!file_out || !(file_out << c))
	throw cado::error("error:: try to write in the file {}.", pathname_out);

    return EXIT_SUCCESS;
}
