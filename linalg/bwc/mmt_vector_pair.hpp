#ifndef MMT_VECTOR_PAIR_HPP_
#define MMT_VECTOR_PAIR_HPP_

#include "matmul_top.h"

class mmt_vector_pair : public std::vector<mmt_vec> {
    matmul_top_data_ptr mmt;
public:
    mmt_vector_pair(matmul_top_data_ptr, int);
    ~mmt_vector_pair();
    mmt_vector_pair(mmt_vector_pair const &) = delete;
    mmt_vector_pair(mmt_vector_pair &&) = delete;
    mmt_vector_pair& operator=(mmt_vector_pair const &) = delete;
    mmt_vector_pair& operator=(mmt_vector_pair &&) = delete;
    mmt_vec * vectors() { return data(); }
};

#endif	/* MMT_VECTOR_PAIR_HPP_ */
