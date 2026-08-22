module;
#include <cstdint>
#include "las-info.hpp"
#include "sieve-methods.hpp"
#include "threadpool.hpp"
#include "las-forwardtypes.hpp"
#include "las-siever-config.hpp"

export module las_choose_sieve_area;

export extern int never_discard;

export template<special_q_data_class T>
bool choose_sieve_area(
        las_info const & las,
        thread_pool & pool,
        special_q_task const & doing,
        siever_config & conf,
        T & Q,
        uint32_t & J);

export template<special_q_data_class T>
bool choose_sieve_area(
        las_info const & las,
        special_q_task const & doing,
        siever_config & conf,
        T & Q,
        uint32_t & J);
