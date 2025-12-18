#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <array>
#include <memory>
#include <utility>
#include <vector>
#include <type_traits>

#include "bucket-push-update.hpp"
#include "bucket.hpp"
#include "fb-types.hpp"
#include "fb.hpp"
#include "las-auxiliary-data.hpp"
#include "las-bkmult.hpp"
#include "las-config.hpp"
#include "las-fill-in-buckets.hpp"
#include "las-process-bucket-region.hpp"
#include "las-qlattice.hpp"
#include "las-report-stats.hpp"
#include "las-siever-config.hpp"
#include "las-smallsieve.hpp"
#include "las-threads-work-data.hpp"
#include "las-where-am-i-proxy.hpp"
#include "las-where-am-i.hpp"
#include "macros.h"
#include "multityped_array.hpp"
#include "tdict.hpp"
#include "threadpool.hpp"
#include "utils_cxx.hpp"
#include "verbose.h"

/* is this in the std library or not ? */
template <typename T> static inline T const & const_ref(T & x)
{
    return x;
}

/***************************************************************************/
/********        Main bucket sieving functions                    **********/

/* {{{ Big question: shall we enable bucket-sieving for powers ?
 *
 * There are several difficulties, in fact. One rationale that yields a
 * straight "no" answer is that such primes make very little difference
 * to the smooth part, so we'd better skip them anyway.
 *
 * But it's not the hardest thing.
 *
 * For the small sieve, we create the small_sieve_data from the factor
 * base entries, and we compute the logp accordingly, per entry.
 *
 * For the bucket-sieve, we use the fact that the factor base is sorted
 * in increasing log(p) order, and we create slices with ranges of primes
 * that have the same round(scale*log(p)).
 *
 * Currently, the factor base is sorted by q=p^k. A power that makes the
 * p-valuation go from p^k0 to p^k1 contributes
 * round(k1*log(p))-round(k0*log(p)). Therefore, sorting by q does not
 * mean that log(p)'s are sorted, and we're in trouble because when we
 * take powers aboard in a slice, their log(p) value is not correctly
 * represented.
 *
 * Previously, we had the behaviour of setting powlim to bucket_thresh-1,
 * effectively preventing powers from appearing in the bucket-sieve.
 *
 * Now powlim is a factor base parameter, and bucket_thresh comes later,
 * so such a default does not work.
 *
 * The strategy we take here is that *if* we see powers here (and we know
 * that will happen only for the fairly rare general entries), then we do
 * something special:
 *  - either we say that we skip over this entry
 *  - or we defer to apply_buckets the expensive computation of a proper
 *    logp value.
 *
 * Currently we do the former. The latter would be viable since only a
 * small fraction of the apply_one_bucket time is devoted to dealing with
 * general entries, so we could imagine having a branch in there for
 * dealing with them. But that would be quite painful. Furthermore, it
 * would then be mandatory to split the entries with same q, same p, but
 * different k0,k1 pairs (we do encounter these), so that the hint would
 * still open the possibility to infer the value of log(p).
 *
 *
 * Note that it would not be possible to sidestep the issue by sorting
 * the vectors of entries by (k1-k0)*log(p) (which would make a
 * difference only for general entries anyway). This is because even
 * sorting by increasing (k1-k0)*log(p) does not guarantee that
 * round(s*k1*log(p))-round(s*k0*log(p)) increases. (counter-example:
 * s=1, k1*log(p)=0.51, k0*log(p)=0.49 diff=0.02 round-round=1
 *      k1*log(p)=1.49, k0*log(p)=0.51 diff=0.98 round-round=0
 * )
 * }}} */
template <class FB_ENTRY_TYPE>
static inline bool discard_power_for_bucket_sieving(FB_ENTRY_TYPE const &)
{
    /* the entry is not a general entry, therefore k is a const thing
     * equal to 1.
     */
    return false;
}
#ifndef BUCKET_SIEVE_POWERS
template <>
inline bool
discard_power_for_bucket_sieving<fb_entry_general>(fb_entry_general const & e)
{
    return e.k > 1;
}
#endif

/* {{{ */
template <int LEVEL> class fill_in_buckets_parameters : public task_parameters
{
  public:
    nfs_work & ws;
    nfs_aux & aux;
    ALGO::special_q_data const & Q;
    int const side;
    fb_slice_interface const * slice;
    plattices_vector_t * plattices_vector; // content changed during fill-in
    plattices_dense_vector_t * plattices_dense_vector; // for sublat
    uint32_t const first_region0_index;
    where_am_I w;

    fill_in_buckets_parameters(nfs_work & _ws, nfs_aux & aux,
                               ALGO::special_q_data const & Q, int const _side,
                               fb_slice_interface const * _slice,
                               plattices_vector_t * _platt,
                               plattices_dense_vector_t * _dplatt,
                               uint32_t const _reg0, where_am_I const & w)
        : ws(_ws)
        , aux(aux)
        , Q(Q)
        , side(_side)
        , slice(_slice)
        , plattices_vector(_platt)
        , plattices_dense_vector(_dplatt)
        , first_region0_index(_reg0)
        , w(w)
    {
    }
};

/* short of a better solution. I know some exist, but it seems way
 * overkill to me.
 *
 * This needs constexpr, though... So maybe I could use a more powerful
 * C++11 trick after all.
 */
#define PREPARE_TEMPLATE_INST_NAMES(F, suffix)                                 \
    template <int> struct CADO_CONCATENATE(F, _name) {                         \
    };                                                                         \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 0);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 1);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 2);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 3);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 4);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 5);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 6);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 7);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 8);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 9)

#define PREPARE_TEMPLATE_INST_NAME(F, suffix, k)                               \
    template <> struct CADO_CONCATENATE(F, _name)<k> {                         \
        static constexpr const char * value = #F "<" #k ">" suffix;            \
    }

/* By tweaking the "" argument, it is possible to have these names
 * embody a suffix like " on side ", so that it's possible tu run
 * parametric timer slots.
 */
PREPARE_TEMPLATE_INST_NAMES(fill_in_buckets_one_slice_internal, "");
PREPARE_TEMPLATE_INST_NAMES(downsort, "");
PREPARE_TEMPLATE_INST_NAMES(downsort_tree, " (dispatcher only)");

#define TEMPLATE_INST_NAME(x, y) CADO_CONCATENATE(x, _name)<y>::value

// For internal levels, the fill-in is not exactly the same as for
// top-level, since the plattices have already been precomputed.
template <int LEVEL, typename TARGET_HINT>
static task_result *
fill_in_buckets_one_slice_internal(worker_thread * worker,
                                   task_parameters * _param, int)
{
    auto * param = static_cast<fill_in_buckets_parameters<LEVEL> *>(_param);

    /* Import some contextual stuff */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(param->aux.th[id]);
    timetree_t & timer(param->aux.get_timer(worker));
    ENTER_THREAD_TIMER(timer);
    nfs_work & ws(param->ws);
    ALGO::special_q_data const & Q(param->Q);
    where_am_I & w(taux.w);
    int const side = param->side;
    nfs_work::side_data & wss(ws.sides[side]);

    MARK_TIMER_FOR_SIDE(timer, side);

    // we're declaring the timer here, but really the work happens below
    // in fill_in_buckets_lowlevel. We happen to have access to
    // param->side here, so we use it to provide a nicer timing report.
    CHILD_TIMER(timer,
                TEMPLATE_INST_NAME(fill_in_buckets_one_slice_internal, LEVEL));

    w = param->w;
    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->plattices_vector->get_index());
    WHERE_AM_I_UPDATE(w, N, param->first_region0_index);

    try {
        /* Get an unused bucket array that we can write to */
        /* clearly, reserve_BA() possibly throws. As it turns out,
         * fill_in_buckets_lowlevel<> does not, at least currently. One
         * could imagine that it could throw, so let's wrap it too.
         */
        auto & BA = wss.reserve_BA<LEVEL, TARGET_HINT>(-1);

        /* Fill the buckets */
        try {
            fill_in_buckets_lowlevel<LEVEL>(BA, ws, Q, *param->plattices_vector,
                                            (param->first_region0_index == 0),
                                            w);
        } catch (buckets_are_full & e) {
            wss.release_BA(BA);
            throw e;
        }
        /* Release bucket array again */
        wss.release_BA(BA);
    } catch (buckets_are_full & e) {
        delete param;
        throw e;
    }
    delete param;
    return new task_result;
}

// At top level.
// We need to interleave the root transforms and the FK walk,
// otherwise, we spend all the time waiting for memory.
// Hence the ugly de-templatization.
// At some point, the code should be re-organized, I'm afraid.
template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static task_result *
fill_in_buckets_toplevel_wrapper(worker_thread * worker MAYBE_UNUSED,
                                 task_parameters * _param, int)
{
    auto * param = static_cast<fill_in_buckets_parameters<LEVEL> *>(_param);

    /* Import some contextual stuff */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(param->aux.th[id]);
    timetree_t & timer(param->aux.get_timer(worker));
    nfs_work & ws(param->ws);
    ALGO::special_q_data const & Q(param->Q);
    where_am_I & w(taux.w);
    int const side = param->side;
    nfs_work::side_data & wss(ws.sides[side]);

    ENTER_THREAD_TIMER(timer);
    MARK_TIMER_FOR_SIDE(timer, side);

#ifndef DISABLE_TIMINGS
    /* This is one of the places where helgrind is likely to complain. We
     * use thread-safe statics. Helgrind can't cope with it,
     * unfortunately. So the error is a false positive.
     *
     * https://sourceforge.net/p/valgrind/mailman/message/32434015/
     */
    timetree_t::accounting_child const local_timer_sentry(timer,
                                                          tdict_slot_for_fibt);
#endif

    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->slice->get_index());
    WHERE_AM_I_UPDATE(w, N, 0);

    try {
        /* Get an unused bucket array that we can write to */
        bucket_array_t<LEVEL, TARGET_HINT> & BA =
            wss.reserve_BA<LEVEL, TARGET_HINT>(-1);
        ASSERT(param->slice);
        auto const * sl = dynamic_cast<fb_slice<FB_ENTRY_TYPE> const *>(param->slice);
        ASSERT_ALWAYS(sl != NULL);
        fill_in_buckets_toplevel<LEVEL, FB_ENTRY_TYPE, TARGET_HINT>(
            BA, ws, *sl, Q, param->plattices_dense_vector, w);
        /* Release bucket array again */
        wss.release_BA(BA);
        delete param;
        return new task_result;
    } catch (buckets_are_full const & e) {
        delete param;
        throw e;
    }
}
/* same for sublat */
template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static task_result *
fill_in_buckets_toplevel_sublat_wrapper(worker_thread * worker,
                                        task_parameters * _param, int)
{
    auto * param = static_cast<fill_in_buckets_parameters<LEVEL> *>(_param);

    /* Import some contextual stuff */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(param->aux.th[id]);
    timetree_t & timer(param->aux.get_timer(worker));
    nfs_work & ws(param->ws);
    ALGO::special_q_data const & Q(param->Q);
    where_am_I & w(taux.w);
    int const side = param->side;
    nfs_work::side_data & wss(ws.sides[side]);

    ENTER_THREAD_TIMER(timer);
    MARK_TIMER_FOR_SIDE(timer, side);

#ifndef DISABLE_TIMINGS
    /* This is one of the places where helgrind is likely to complain. We
     * use thread-safe statics. Helgrind can't cope with it,
     * unfortunately. So the error is a false positive.
     *
     * https://sourceforge.net/p/valgrind/mailman/message/32434015/
     */
    timetree_t::accounting_child const local_timer_sentry(timer,
                                                          tdict_slot_for_fibt);
#endif

    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->slice->get_index());
    WHERE_AM_I_UPDATE(w, N, 0);

    try {
        /* Get an unused bucket array that we can write to */
        bucket_array_t<LEVEL, TARGET_HINT> & BA =
            wss.reserve_BA<LEVEL, TARGET_HINT>(-1);
        ASSERT(param->slice);
        fill_in_buckets_toplevel_sublat<LEVEL, FB_ENTRY_TYPE>(
            BA, ws, Q,
            *dynamic_cast<fb_slice<FB_ENTRY_TYPE> const *>(param->slice),
            param->plattices_dense_vector, w);
        /* Release bucket array again */
        wss.release_BA(BA);
        delete param;
        return new task_result;
    } catch (buckets_are_full const & e) {
        delete param;
        throw e;
    }
}

/* Whether or not we want fill_in_buckets_one_slice to be templatized
 * both for LEVEL and n is not clear. At some point, we're doing code
 * bloat for almost nothing.
 *
 * Now given the code below, it's easy enough to arrange so that we go
 * back to the virtual base fb_slice_interface.
 */
template <int LEVEL, typename TARGET_HINT> struct push_slice_to_task_list {
    thread_pool & pool;
    fill_in_buckets_parameters<LEVEL> model;
    push_slice_to_task_list(thread_pool & pool,
                            fill_in_buckets_parameters<LEVEL> const & m)
        : pool(pool)
        , model(m)
    {
    }
    size_t pushed = 0;
    template <typename T> void operator()(T const & s)
    {
        auto * param = new fill_in_buckets_parameters<LEVEL>(model);
        param->slice = &s;
        using entry_t = typename T::entry_t;
        task_function_t f =
            fill_in_buckets_toplevel_wrapper<LEVEL, entry_t, TARGET_HINT>;
        pool.add_task(f, param, 0, 0, s.get_weight());
        pushed++;
    }
};
template <int LEVEL, typename TARGET_HINT>
struct push_slice_to_task_list_saving_precomp {
    thread_pool & pool;
    fb_factorbase::slicing::part const & P;
    fill_in_buckets_parameters<LEVEL> model;
    /* precomp_plattice_dense_t == std::vector<plattices_dense_vector_t> */
    typename precomp_plattice_dense_t<LEVEL>::type & Vpre;
    size_t pushed = 0;
    push_slice_to_task_list_saving_precomp(
        thread_pool & pool, fb_factorbase::slicing::part const & P,
        fill_in_buckets_parameters<LEVEL> const & m,
        typename precomp_plattice_dense_t<LEVEL>::type & Vpre)
        : pool(pool)
        , P(P)
        , model(m)
        , Vpre(Vpre)
    {
    }
    template <typename T> void operator()(T const & s)
    {
        /* we're pushing the global index, relative to all fb parts */
        slice_index_t const idx = s.get_index();
        ASSERT_ALWAYS((size_t)idx == pushed);

        plattices_dense_vector_t & pre(Vpre[idx]);

        auto * param = new fill_in_buckets_parameters<LEVEL>(model);
        param->slice = &s;
        param->plattices_dense_vector = &pre;

        using entry_t = typename T::entry_t;
        task_function_t f =
            fill_in_buckets_toplevel_sublat_wrapper<LEVEL, entry_t,
                                                    TARGET_HINT>;
        pool.add_task(f, param, 0, 0, s.get_weight());
        pushed++;
    }
};

template <int LEVEL, typename TARGET_HINT>
static void fill_in_buckets_one_side(nfs_work & ws, nfs_aux & aux,
                                     ALGO::special_q_data const & Q,
                                     thread_pool & pool, int const side,
                                     where_am_I & w)
{
    timetree_t & timer(aux.rt.timer);
    nfs_work::side_data & wss(ws.sides[side]);

    /* We're just pushing tasks, here. */
    BOOKKEEPING_TIMER(timer);

    fill_in_buckets_parameters<LEVEL> const model(ws, aux, Q, side, NULL, NULL,
                                                  NULL, 0, w);

    auto const & BA_ins = wss.bucket_arrays<LEVEL, TARGET_HINT>();

    verbose_fmt_print(0, 3,
            "# Filling the side-{} {}{} buckets ({} groups of {} buckets)\n",
            side,
            LEVEL, TARGET_HINT::rtti[0],
            BA_ins.size(), BA_ins[0].n_bucket);

    /* We'd like to also display info on the slices we're about to run
     * FIB on, but the multityped_array interface won't let us do this
     * easily.
     */

    if (!Q.sublat.m) {
        /* This creates a task meant to call
         * fill_in_buckets_toplevel_wrapper */
        push_slice_to_task_list<LEVEL, TARGET_HINT> F(pool, model);
        wss.fbs->get_part(LEVEL).foreach_slice(F);
    } else {
        /* This creates a task meant to call
         * fill_in_buckets_toplevel_sublat_wrapper */
        auto & Vpre(wss.precomp_plattice_dense.get<LEVEL>());
        fb_factorbase::slicing::part const & P = wss.fbs->get_part(LEVEL);

        /* This way we can spare the need to expose the copy contructor
         * of the container's value_type */
        if (Q.sublat.i0 == 0 && Q.sublat.j0 == 1) {
            /* first sublat */
            Vpre = typename precomp_plattice_dense_t<LEVEL>::type(P.nslices());
        }

        ASSERT_ALWAYS(Vpre.size() == P.nslices());
        push_slice_to_task_list_saving_precomp<LEVEL, TARGET_HINT> F(
            pool, P, model, Vpre);
        P.foreach_slice(F);
    }
}


/* This is a compile-time loop over the possible values from 1 to level,
 * and 0 errors out. */
template<int level, typename hint_t>
struct fib1s_caller_s : public fib1s_caller_s<level-1, hint_t> {
    template<typename... Args>
    void operator()(nfs_work & ws, Args&& ...args) const {
        if (ws.toplevel == level)
            fill_in_buckets_one_side<level, hint_t>(ws, std::forward<Args>(args)...);
        else
            fib1s_caller_s<level-1, hint_t>::operator()(ws, std::forward<Args>(args)...);
    }
};
template<typename hint_t>
struct fib1s_caller_s<0, hint_t> {
    template<typename... Args>
    void operator()(nfs_work &, Args&& ...) const {
        ASSERT_ALWAYS(0);
    }
};

template<int level, typename hint_t, typename... Args>
inline void fib_one_side(nfs_work & ws, Args&& ...args)
{
    fib1s_caller_s<level, hint_t>()(ws, std::forward<Args>(args)...);
}

void fill_in_buckets_toplevel_multiplex(nfs_work & ws, nfs_aux & aux,
        ALGO::special_q_data const & Q, thread_pool & pool, int side, where_am_I & w)
{
    // per se, we're not doing anything here.
    // CHILD_TIMER(timer, __func__);
    if (ws.conf.needs_resieving()) {
        fib_one_side<MAX_TOPLEVEL, shorthint_t>(ws, aux, Q, pool, side, w);
    } else {
        fib_one_side<MAX_TOPLEVEL, emptyhint_t>(ws, aux, Q, pool, side, w);
    }
}

/* }}} */

#ifdef SIQS_SIEVE
#include "siqs-fill-in-buckets.inl"
#else
#include "las-fill-in-buckets.inl"
#endif

/* Instances to be compiled */

// A fake level 0, to avoid infinite loop during compilation.
template <>
void downsort_tree<0>(nfs_work &,
                      std::shared_ptr<nfs_work_cofac>,
                      std::shared_ptr<nfs_aux>,
                      ALGO::special_q_data const &,
                      thread_pool &, uint32_t,
                      uint32_t,
                      std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> &,
                      where_am_I &)
{
    ASSERT_ALWAYS(0);
}

// Now the exported instances

template void downsort_tree<1>(
    nfs_work &, std::shared_ptr<nfs_work_cofac>, std::shared_ptr<nfs_aux> aux_p,
    ALGO::special_q_data const &, thread_pool &, uint32_t, uint32_t,
    std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> &, where_am_I &);


#if MAX_TOPLEVEL >= 3
template void
downsort_tree<2>(nfs_work &, std::shared_ptr<nfs_work_cofac>,
                 std::shared_ptr<nfs_aux>, ALGO::special_q_data const & Q,
                 thread_pool &, uint32_t, uint32_t,
                 std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> &,
                 where_am_I &);
#endif
static_assert(MAX_TOPLEVEL == 3);
