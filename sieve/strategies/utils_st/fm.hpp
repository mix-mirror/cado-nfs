#ifndef CADO_FM_HPP
#define CADO_FM_HPP

#include <istream>
#include <ostream>
#include <vector>

#include "fmt/base.h"
#include "fmt/ostream.h"

#include "facul_ecm.h"
#include "facul_method.hpp"

/* A factoring method, together with what a benchmark measured about it.
 *
 * Only the method family and the bounds B1 and B2 are part of the
 * measured data. The ECM parameter (sigma) is not: it is drawn when the
 * method is instantiated for a bench or for a real factorization, so
 * params.parameter is left at 0 here and carries no meaning.
 */
struct factoring_method {
    facul_method::parameters params;

    /* proba[i] is the probability that this method finds a prime factor
     * of len_p_min + i bits. */
    unsigned int len_p_min = 0;
    std::vector<double> proba;

    /* time[i] is the time this method takes on an input whose size puts
     * it in bucket i; see time_index() below. */
    std::vector<double> time;

    factoring_method() = default;
    explicit factoring_method(facul_method::parameters p)
        : params(std::move(p))
    {
    }

    /* Build the parameters from the four numbers that the file format
     * stores, which is also how the generator passes them around. */
    static facul_method::parameters
    make_parameters(facul_method_code method, ec_parameterization_t curve,
                    unsigned long B1, unsigned long B2);

    /* from the four numbers of one record of the file format */
    static factoring_method from_fields(unsigned long method,
                                        unsigned long curve, unsigned long B1,
                                        unsigned long B2)
    {
        return factoring_method(make_parameters(facul_method_code(method),
                                                ec_parameterization_t(curve),
                                                B1, B2));
    }

    facul_method_code method() const { return params.method; }
    unsigned long B1() const { return params.B1; }
    unsigned long B2() const { return params.B2; }

    /* 0 for everything but ECM, matching what the file format stores. */
    ec_parameterization_t curve() const
    {
        return params.method == EC_METHOD ? params.parameterization
                                          : ec_parameterization_t(0);
    }

    /* The generator uses methods with B1 == B2 == 0 as placeholders for
     * "no method in this slot". */
    bool is_zero() const { return params.is_null(); }
    void put_zero();

    /* Only the method and its bounds are compared, not the measurements:
     * this answers "are these the same method?". */
    bool same_method_as(factoring_method const & o) const;

    /* time[] is indexed by the number of half-words of the input, the
     * way bench_time() measures it: bucket 0 covers anything up to one
     * word, then one bucket per extra half-word. Returns an index that
     * may be past the end of a short time[]; callers clamp. */
    static unsigned int time_index(unsigned int r);

    /* time[] entry to use for an r-bit input, clamped to what we have. */
    double time_for(unsigned int r) const;
};

std::istream & operator>>(std::istream & is, factoring_method &);
std::ostream & operator<<(std::ostream & os, factoring_method const &);

namespace fmt
{
template <> struct formatter<factoring_method> : ostream_formatter {
};
} // namespace fmt

#endif /* CADO_FM_HPP */
