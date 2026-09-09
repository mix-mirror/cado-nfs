#include "cado.h" // IWYU pragma: keep

#include "fmt/base.h"

#include "cxx_mpz.hpp"
#include "getprime.h"
#include "params.hpp"

struct prog
{
    parameter_mandatory<cxx_mpz, "N",
        "the integer to check"> n;
    parameter_mandatory<unsigned long, "B",
        "check divisibility by square of prime smaller than this bound"> B;
    parameter_switch<"skip2", "do not check divisibility by 2^2"> skip2;

    static void configure(cxx_param_list & pl)
    {
        decltype(prog::n)::configure(pl);
        decltype(prog::B)::configure(pl);
    }

    explicit prog(cxx_param_list & pl)
        : n(pl), B(pl)
    {
    }

    bool is_squarefree_up_to() const
    {
        cxx_mpz c(n);
        mpz_abs(c, c);
        for (auto const p: prime_range(skip2 ? 3u : 2u, B)) {
            if (mpz_divisible_ui_p(c, p)) {
                mpz_divexact_ui(c, c, p);
                if (mpz_divisible_ui_p(c, p)) {
                    fmt::print(stderr,
                               "{} has a square factor {}^2 below {}\n",
                               n.value, p, B.value);
                    return false;
                }
                if (c == 1u)
                    break;
            }
        }
        return true;
    }
};

int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    prog::configure(pl);
    pl.process_command_line(argc, argv, true);

    prog P(pl);

    pl.print_command_line(stdout);
    fflush(stdout);

    if (pl.warn_unused())
        pl.fail("Error, unused parameters are given\n");

    return P.is_squarefree_up_to() ? EXIT_SUCCESS : EXIT_FAILURE;
}
