#ifndef MMT_VECTOR_PAIR_HPP_
#define MMT_VECTOR_PAIR_HPP_

#include "matmul_top_vec.hpp"
#include <tuple>

struct matmul_top_matrix;
struct matmul_top_data;

#include "matmul_top.hpp"

class mmt_vector_pair : private std::vector<mmt_vec> {
    /* if direction = 0 we have:
     * [1] = [0] * M0
     * [2] = [1] * M1
     * [...]
     * [n] = [n-1] * M_{n-1} (if n is odd)
     * [0] = [n-1] * M_{n-1} (if n is even)
     *
     * if direction = 1 we have:
     * [1] = M_{n-1} * [0]
     * [2] = M_{n-2} * [1]
     * [...]
     * [n] = M_0 * [n-1] (if n is odd)
     * [0] = M_0 * [n-1] (if n is even)
     *
     */
    bool nmats_odd;
public:
    int direction;
    inline mmt_vec const & input_vector() const { return front(); }
    inline mmt_vec & input_vector() { return front(); }
    inline mmt_vec const & output_vector() const { return nmats_odd ? back() : front(); }
    inline mmt_vec & output_vector() { return nmats_odd ? back() : front(); }
    inline mmt_vec & operator[](int i) { return ((std::vector<mmt_vec> &)*this)[i]; }
    inline mmt_vec const & operator[](int i) const { return ((std::vector<mmt_vec> const &)*this)[i]; }

    mmt_vector_pair(matmul_top_data &, int);
    mmt_vector_pair(mmt_vector_pair const &) = delete;
    mmt_vector_pair(mmt_vector_pair &&) = delete;
    mmt_vector_pair& operator=(mmt_vector_pair const &) = delete;
    mmt_vector_pair& operator=(mmt_vector_pair &&) = delete;
    // mmt_vec * vectors() { return data(); }
    // mmt_vec * vectors() { return reinterpret_cast<mmt_vec*>(data()); }

    struct multiplication_chain {
        matmul_top_data const & mmt;
        mmt_vector_pair & vecs;
        multiplication_chain(matmul_top_data const & mmt, mmt_vector_pair & vecs)
            : mmt(mmt)
            , vecs(vecs)
        {
        }
        typedef std::tuple<
                matmul_top_matrix const &, /* matrix */
                mmt_vec &,                 /* input */
                mmt_vec &,                 /* output */
                bool,                      /* whether output is vec[n] */
                size_t>                    /* matrix index*/
                value_type_base;

        class value_type : private value_type_base
        {
            value_type_base& super() { return *this; }
            public:
            template<typename... Args> value_type(Args&&... args)
                : value_type_base(std::forward<Args>(args)...) {}
            matmul_top_matrix const & matrix() { return std::get<0>(super()); }
            mmt_vec & src() { return std::get<1>(super()); }
            mmt_vec & dst() { return std::get<2>(super()); }
            bool output_is_extra_vector() { return std::get<3>(super()); }
            int matrix_index() { return std::get<4>(super()); }
        };

        value_type get(size_t i) {
            return value_type(
                    mmt.matrices[vecs.direction ? mmt.matrices.size()-1-i : i],
                    vecs[i],
                    vecs[(i == mmt.matrices.size() - 1 && !vecs.nmats_odd)
                    ? 0 : (i + 1)],
                    i == mmt.matrices.size() - 1 && vecs.nmats_odd,
                    vecs.direction ? mmt.matrices.size()-1-i : i
                    );
        }
        struct sentinel {
            multiplication_chain & c;
            size_t i;
            sentinel(multiplication_chain & c, size_t i) : c(c), i(i) {}
            inline sentinel& operator++() { ++i; return *this; }
            inline bool operator!=(sentinel const & e) { return i != e.i; }
            multiplication_chain::value_type operator*() const {
                return c.get(i);
            }
        };
        sentinel begin() { return sentinel(*this, 0); }
        sentinel end() { return sentinel(*this, mmt.matrices.size()); }
    };

    friend struct multiplication_chain;
    multiplication_chain chain(matmul_top_data & mmt) {
        return multiplication_chain(mmt, *this);
    }
};

#endif	/* MMT_VECTOR_PAIR_HPP_ */
