#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <sstream>

#include "fmt/base.h"

#include "decomp.hpp"
#include "macros.h"
#include "tab_decomp.hpp"

int main()
{
    decomp el1 {1000, { 1U, 2U, 3U} };
    tabular_decomp t {el1};
    if (!(el1 == t[0])) {
        fmt::print(stderr, "error with the test(1)!!!\n");
        return EXIT_FAILURE;
    }
    t.push_back(el1);
    // test realloc()
    const decomp el2 { 10000, { 11U, 10U, 9U, 8U, 7U }};
    t.push_back(el2);
    if (t[2] != el2) {
        fmt::print(stderr, "error with the test(2)!!!\n");
        return EXIT_FAILURE;
    }
    if (t[2] == el1) {
        fmt::print(stderr, "error with the test(3)!!!\n");
        return EXIT_FAILURE;
    }
    // set and get
    el1 = el2;
    if (el1 != el2) {
        fmt::print(stderr, "error with the test(4)!!!\n");
        return EXIT_FAILURE;
    }
    // round trip through the stream operators
    {
        std::ostringstream os;
        os << t;

        tabular_decomp t2;
        std::istringstream is(os.str());
        if (!(is >> t2)) {
            fmt::print(stderr, "read error on what we just wrote\n");
            return EXIT_FAILURE;
        }
        if (t2.size() != t.size()) {
            fmt::print(stderr, "error with the test(5)!!!\n");
            return EXIT_FAILURE;
        }
        for (size_t i = 0; i < t.size(); i++) {
            if (t2[i] != t[i]) {
                fmt::print(stderr, "error with the test(5)!!!\n");
                return EXIT_FAILURE;
            }
        }
        /* nb_elem is part of what a decomp compares as */
        t2[1].nb_elem = 21;
        if (t2[1] == t[1]) {
            fmt::print(stderr, "error with the test(6)!!!\n");
            return EXIT_FAILURE;
        }
    }

    /* a malformed record must fail the stream, not abort */
    {
        std::istringstream is("[ 1 2 3 ; 1000\n");
        tabular_decomp t3;
        if (is >> t3) {
            fmt::print(stderr, "error with the test(7)!!!\n");
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}
