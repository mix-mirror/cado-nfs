#include "cado.h" // IWYU pragma: keep

#include <cmath>

#include <istream>
#include <ostream>
#include <utility>

#include "arith/modredc_ul.h" // MODREDCUL_MAXBITS

#include "fm.hpp"
#include "utils_cxx.hpp"

facul_method::parameters
factoring_method::make_parameters(facul_method_code method,
                                  ec_parameterization_t curve,
                                  unsigned long B1, unsigned long B2)
{
    /* The file format has a slot for the curve next to every method,
     * even those that have no curve. gfm fills it with whichever
     * parameterization its sweep was on, so it is meaningless for PM1
     * and PP1 -- but it is carried through verbatim, because that is
     * what the format stores and what the readers of these files have
     * always seen. */
    if (method == EC_METHOD) {
        return {method, B1, B2, curve, 0, 1};
    }
    facul_method::parameters p {method, B1, B2};
    p.parameterization = curve;
    return p;
}

void factoring_method::put_zero()
{
    params.B1 = 0;
    params.B2 = 0;
    for (auto & x: proba)
        x = 0;
    for (auto & x: time)
        x = 0;
}

bool factoring_method::same_method_as(factoring_method const & o) const
{
    return params.method == o.params.method &&
           params.parameterization == o.params.parameterization &&
           params.B1 == o.params.B1 && params.B2 == o.params.B2;
}

unsigned int factoring_method::time_index(unsigned int r)
{
    /* We add 0.5 to the length of one word, because our times are
     * measured for an inclusive length: with MODREDCUL_MAXBITS = 64, a
     * cofactor fits in one word when its length is less than *or equal
     * to* 64 bits. Without the 0.5 the equality case would land in the
     * next bucket. */
    double const half_word = (MODREDCUL_MAXBITS + 0.5) / 2.0;
    unsigned int const number_half_wd = (unsigned int) floor(r / half_word);
    return number_half_wd < 2 ? 0 : number_half_wd - 1;
}

double factoring_method::time_for(unsigned int r) const
{
    unsigned int const i = time_index(r);
    return i >= time.size() ? time.back() : time[i];
}

std::ostream & operator<<(std::ostream & os, factoring_method const & fm)
{
    /* This has to stay byte-compatible with what gst and benchfm read:
     * four numbers, then the probabilities behind the bit size they
     * start at, then the timings, each field closed by a '|'. */
    os << fm.params.method << ' ' << fm.params.parameterization << ' '
       << fm.params.B1 << ' ' << fm.params.B2 << " | ";
    os << fm.len_p_min << ' ';
    for (double const p: fm.proba)
        os << fmt::format("{:f} ", p);
    os << "| ";
    for (double const t: fm.time)
        os << fmt::format("{:f} ", t);
    os << "|\n";
    return os;
}

/* read doubles until the next '|' */
static bool read_doubles(std::istream & is, std::vector<double> & v)
{
    for (;;) {
        is >> std::ws;
        if (!is.good() || is.peek() == '|')
            break;
        double x;
        if (!(is >> x))
            return false;
        v.push_back(x);
    }
    return true;
}

std::istream & operator>>(std::istream & is, factoring_method & fm)
{
    unsigned long m[4];
    for (auto & x: m)
        if (!(is >> x))
            return is;

    if (!(is >> std::ws >> expect("|")))
        return is;

    unsigned int len_p_min;
    if (!(is >> len_p_min))
        return is;

    std::vector<double> proba;
    std::vector<double> time;

    if (!read_doubles(is, proba) || !(is >> std::ws >> expect("|")))
        return is;
    if (!read_doubles(is, time) || !(is >> std::ws >> expect("|")))
        return is;

    fm.params = factoring_method::make_parameters(
        facul_method_code(m[0]), ec_parameterization_t(m[1]), m[2], m[3]);
    fm.len_p_min = len_p_min;
    fm.proba = std::move(proba);
    fm.time = std::move(time);
    return is;
}
