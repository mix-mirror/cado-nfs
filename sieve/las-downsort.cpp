#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "bucket.hpp"
#include "fb-types.hpp"
#include "las-auxiliary-data.hpp"
#include "las-bkmult.hpp"
#include "las-config.hpp"
#include "las-downsort.hpp"
#include "las-fill-in-buckets.hpp"
#include "las-globals.hpp"
#include "las-process-bucket-region.hpp"
#include "las-report-stats.hpp"
#include "las-siever-config.hpp"
#include "las-threads-work-data.hpp"
#include "las-where-am-i-proxy.hpp"
#include "macros.h"
#include "multityped_array.hpp"
#include "las-plattice.hpp"
#include "threadpool.hpp"
#include "las-where-am-i-proxy.hpp"
#include "las-where-am-i.hpp"

#ifdef SIQS_SIEVE
#include "siqs-fill-in-buckets.inl"
#else
#include "las-fill-in-buckets.inl"
#endif

/* multithreaded implementation of the downsort procedure. It becomes a
 * bottleneck sooner than one might think.
 *
 */

/* This is auxiliary only. We downsort stuff that we already downsorted.
 * So it applies only if LEVEL+1 is itself not the toplevel.
 * For this reason, we must have a specific instantiation that reduces
 * this to a no-op if LEVEL+1>=3, because there's no longhint_t for level
 * 3 presently.
 */

template <int LEVEL, bool WITH_HINTS>
static void downsort_aux(nfs_work & ws,
                         nfs_aux & aux, thread_pool & pool, task_group & tg, int side,
                         uint32_t bucket_index, where_am_I & w)
{
    static_assert(LEVEL <= MAX_TOPLEVEL - 1);

    using my_longhint_t = hints_proxy<WITH_HINTS>::l;

    nfs_work::side_data & wss(ws.sides[side]);


    auto const & BA_ins = wss.bucket_arrays<LEVEL + 1, my_longhint_t>();
    auto & BA_outs = wss.bucket_arrays<LEVEL, my_longhint_t>();
    ASSERT_ALWAYS(BA_ins.size() == BA_outs.size());

    verbose_fmt_print(0, 3,
            "# Downsorting the side-{} {}{} buckets ({} groups of {} buckets"
            ", taking only bucket {}/{})"
            " to {}{} buckets ({} groups of {} buckets)\n",
            side,
            LEVEL + 1, my_longhint_t::rtti[0],
            BA_ins.size(), BA_ins[0].n_bucket,
            bucket_index, BA_ins[0].n_bucket,
            LEVEL, my_longhint_t::rtti[0],
            BA_outs.size(), BA_outs[0].n_bucket);


    // What comes from already downsorted data above:
    for (auto const & BA_in: BA_ins) {
        pool.add_task(
            tg,
            [&aux, &ws, side, bucket_index, &BA_in](worker_thread * worker, where_am_I && w) {
                nfs_work::side_data & wss(ws.sides[side]);
                fb_factorbase::slicing const & fbs(*wss.fbs);
                int const id = worker->rank();
                nfs_aux::thread_data & taux(aux.th[id]);
                timetree_t & timer(aux.get_timer(worker));
                ENTER_THREAD_TIMER(timer);
                MARK_TIMER_FOR_SIDE(timer, side);
                taux.w = std::move(w);
                CHILD_TIMER(timer, fmt::format("downsort<{}>", LEVEL));
                auto tt = worker->trace(chronograms::DS(side,LEVEL,bucket_index));
                downsort<LEVEL + 1>(fbs,
                        wss.acquire_BA<LEVEL, my_longhint_t>(wss.rank_BA(BA_in)),
                        BA_in, bucket_index, taux.w);
            }, where_am_I(w));
    }
}

#if MAX_TOPLEVEL == 2
template <>
void downsort_aux<1, false>(nfs_work &, nfs_aux &,
                     thread_pool &, task_group &, int, uint32_t, where_am_I &)
{
}
template <>
void downsort_aux<1, true>(nfs_work &, nfs_aux &,
                     thread_pool &, task_group &, int, uint32_t, where_am_I &)
{
}
#endif
#if MAX_TOPLEVEL == 3
template <>
void downsort_aux<2, false>(nfs_work &, nfs_aux &,
                     thread_pool &, task_group &, int, uint32_t, where_am_I &)
{
}
template <>
void downsort_aux<2, true>(nfs_work &, nfs_aux &,
                     thread_pool &, task_group &, int, uint32_t, where_am_I &)
{
}
#endif
static_assert(MAX_TOPLEVEL == 3);

// For internal levels, the fill-in is not exactly the same as for
// top-level, since the plattices have already been precomputed.
template <int LEVEL, typename TARGET_HINT>
static void
fill_in_buckets_one_slice_internal(worker_thread * worker,
        nfs_work & ws,
        nfs_aux & aux,
        ALGO::special_q_data const & Q,
        int side,
        fb_slice_interface *,
        plattices_vector_t * plattices_vector,
        plattices_dense_vector_t *,
        uint32_t const first_region0_index)
{
    static_assert(!TARGET_HINT::is_long_v);

    /* Import some contextual stuff */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(aux.th[id]);
    timetree_t & timer(aux.get_timer(worker));
    ENTER_THREAD_TIMER(timer);
    where_am_I & w(taux.w);
    nfs_work::side_data & wss(ws.sides[side]);

    MARK_TIMER_FOR_SIDE(timer, side);

    // we're declaring the timer here, but really the work happens below
    // in fill_in_buckets_lowlevel. We happen to have access to
    // param->side here, so we use it to provide a nicer timing report.
    CHILD_TIMER(timer,
            fmt::format("fill_in_buckets_one_slice_internal<{}>", LEVEL));

    WHERE_AM_I_UPDATE(w, side, side);
    WHERE_AM_I_UPDATE(w, i, plattices_vector->get_index());
    WHERE_AM_I_UPDATE(w, N, first_region0_index);

    try {
        auto acquired = wss.reserve_BA<LEVEL, TARGET_HINT>();
        auto tt = worker->trace(chronograms::FIB(
                    side,
                    LEVEL,
                    wss.rank_BA(acquired.access()),
                    plattices_vector->get_index()));
            
        /* Get an unused bucket array that we can write to */
        /* clearly, reserve_BA() possibly throws. As it turns out,
         * fill_in_buckets_lowlevel<> does not, at least currently. One
         * could imagine that it could throw, so let's wrap it too.
         */
        fill_in_buckets_lowlevel<LEVEL, TARGET_HINT>(
                acquired.access(),
                ws, Q, *plattices_vector,
                first_region0_index, w);
    } catch (buckets_are_full & e) {
        e.side = side;
        throw e;
    }
}

// first_region0_index is a way to remember where we are in the tree.
// The depth-first is a way to process all the the regions of level 0 in
// increasing order of j-value.
// first_region0_index * nb_lines_per_region0 therefore gives the j-line
// where we are. This is what is called N by WHERE_AM_I and friends.

template <int LEVEL, bool WITH_HINTS>
static void downsort_tree_inner(
    nfs_work & ws,
    std::shared_ptr<nfs_work_cofac> wc_p,
    std::shared_ptr<nfs_aux> aux_p,
    ALGO::special_q_data const & Q,
    thread_pool & pool,
    uint32_t bucket_index, /* for the current level ! */
    uint32_t first_region0_index,
    std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> & precomp_plattices,
    where_am_I & w)
{
    /* LEVEL is not the toplevel here, so we must have the following: */
    static_assert(LEVEL <= MAX_TOPLEVEL - 1);

    int const nsides = ws.las.cpoly.nsides();
    nfs_aux & aux(*aux_p);
    timetree_t & timer(aux.rt.timer);

    using my_longhint_t = hints_proxy<WITH_HINTS>::l;
    using my_shorthint_t = hints_proxy<WITH_HINTS>::s;

    CHILD_TIMER(timer, fmt::format("downsort_tree<{}>", LEVEL));
    TIMER_CATEGORY(timer, sieving_mixed());
    ASSERT_ALWAYS(LEVEL > 0);

    WHERE_AM_I_UPDATE(w, N, first_region0_index);

    std::vector<task_group> ds_tgs(nsides);
    std::vector<task_group> ds_aux_tgs(nsides);
    std::vector<task_group> fib_tgs(nsides);
    std::vector<task_group> sss_tgs(nsides);

    for (int side = 0; side < nsides; ++side) {
        nfs_work::side_data & wss(ws.sides[side]);
        if (wss.no_fb())
            continue;

        WHERE_AM_I_UPDATE(w, side, side);
        TIMER_CATEGORY(timer, sieving(side));

        auto & ds_tg = ds_tgs[side];
        auto & ds_aux_tg = ds_aux_tgs[side];
        auto & fib_tg = fib_tgs[side];

        /* FIRST: Downsort what is coming from the level above, for this
         * bucket index */
        // All these BA are global stuff; see reservation_group.
        // We reserve those where we write, and access the ones for
        // reading without reserving. We require that things at level
        // above are finished before entering here.

        {
            auto const & BA_ins = wss.bucket_arrays<LEVEL + 1, my_shorthint_t>();
            auto & BA_outs = wss.bucket_arrays<LEVEL, my_longhint_t>();
            /* otherwise the code here can't work */
            /* see also the comment above in downsort_aux */
            ASSERT_ALWAYS(BA_ins.size() == BA_outs.size());

            verbose_fmt_print(0, 3,
                    "# Downsorting the side-{} {}{} buckets ({} groups of {} buckets"
                    ", taking only bucket {}/{})"
                    " to {}{} buckets ({} groups of {} buckets)\n",
                    side,
                    LEVEL + 1, my_shorthint_t::rtti[0],
                    BA_ins.size(), BA_ins[0].n_bucket,
                    bucket_index, BA_ins[0].n_bucket,
                    LEVEL, my_longhint_t::rtti[0],
                    BA_outs.size(), BA_outs[0].n_bucket);

            ASSERT_ALWAYS(BA_ins.size() == BA_outs.size());

            wss.reset_all_pointers<LEVEL, my_longhint_t>();

            for (auto const & BA_in: BA_ins) {
                pool.add_task(
                    ds_tg,
                    [&aux, &ws, side, &BA_in, bucket_index](worker_thread * worker, where_am_I && w) {
                        nfs_work::side_data & wss(ws.sides[side]);
                        fb_factorbase::slicing const & fbs(*wss.fbs);
                        int const id = worker->rank();
                        nfs_aux::thread_data & taux(aux.th[id]);
                        timetree_t & timer(aux.get_timer(worker));
                        taux.w = std::move(w);
                        ENTER_THREAD_TIMER(timer);
                        MARK_TIMER_FOR_SIDE(timer, side);
                        CHILD_TIMER(timer, fmt::format("downsort<{}>", LEVEL));
                        auto tt = worker->trace(chronograms::DS(side, LEVEL, bucket_index));

                        downsort<LEVEL + 1>(fbs,
                                wss.acquire_BA<LEVEL, my_longhint_t>(wss.rank_BA(BA_in)),
                                BA_in, bucket_index,
                                taux.w);
                    },
                    where_am_I(w));
            }
        }

        /* Once ds_tg completes for this side, enqueue all FIB slice tasks into fib_tg */
        ds_tg.on_complete(
            [&ws, &aux, &pool, &precomp_plattices, &ds_aux_tg, &fib_tg, side, first_region0_index, &Q, bucket_index, w_val = where_am_I(w)]() mutable {
                /* SECOND: fill in buckets at this level, for this region. */
                /* We do so only once ds_tg completes for this side */
                auto do_fib = [&ws, &aux, &pool, &precomp_plattices, &fib_tg, side, first_region0_index, &Q]() {
                    nfs_work::side_data & wss(ws.sides[side]);
                    wss.reset_all_pointers<LEVEL, my_shorthint_t>();

                    auto & BA_outs = wss.bucket_arrays<LEVEL, my_shorthint_t>();
                    auto & lattices = precomp_plattices[side].get<LEVEL>();

                    verbose_fmt_print(0, 3,
                            "# Filling the side-{} {}{} buckets ({} groups of {} buckets)"
                            " using {} precomputed lattices\n",
                            side,
                            LEVEL, my_shorthint_t::rtti[0],
                            BA_outs.size(), BA_outs[0].n_bucket,
                            lattices.size());
                    if (!lattices.empty()) {
                        verbose_fmt_print(0, 3,
                                "#   lattices go from slice {} ({} primes) to slice {} ({} primes)\n",
                                lattices.front().get_index(), lattices.front().size(),
                                lattices.back().get_index(), lattices.back().size()
                            );
                    }

                    for (auto & it: lattices) {
                        pool.add_task(fib_tg, thread_pool::QUEUE_GENERIC, it.get_weight(),
                                fill_in_buckets_one_slice_internal<LEVEL, my_shorthint_t>,
                                std::ref(ws), std::ref(aux), std::ref(Q), side, nullptr, &it, nullptr, first_region0_index);
                    }
                };

                if (LEVEL < ws.toplevel - 1) {
                    downsort_aux<LEVEL, WITH_HINTS>(ws, aux, pool, ds_aux_tg, side, bucket_index, w_val);
                    ds_aux_tg.on_complete(do_fib);
                } else {
                    do_fib();
                }
            });
    }

    if (LEVEL == 1) {
        /* Prepare for PBR: we need to precompute the small sieve positions
         * for all the small sieved primes.
         *
         * For ws.toplevel==1, we don't reach here, of course, and the
         * corresponding initialization is done with identical code in
         * las.cpp
         */
        ASSERT(ws.toplevel > 1);
        for (int side = 0; side < nsides; side++) {
            nfs_work::side_data const & wss(ws.sides[side]);
            if (wss.no_fb())
                continue;

            auto & sss_tg = sss_tgs[side];
            wss.ssd->small_sieve_prepare_many_start_positions(
                    pool, &sss_tg,
                    first_region0_index,
                    std::min(SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE, ws.nb_buckets[1]),
                    ws.conf.logI, Q.sublat);
            sss_tg.on_complete([&wss]() {
                wss.ssd->small_sieve_activate_many_start_positions();
            });
        }
    }

    for (int side = 0; side < nsides; side++) {
        nfs_work::side_data const & wss(ws.sides[side]);
        if (wss.no_fb())
            continue;
        ds_tgs[side].wait();
        ds_aux_tgs[side].wait();
        fib_tgs[side].wait();
        if (LEVEL == 1)
            sss_tgs[side].wait();
    }

    /* RECURSE */
    if (LEVEL > 1) {
        size_t const(&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
        verbose_fmt_print(0, 3,
                "# recursively downsort level-{} buckets ({} buckets)"
                " to level {} (+ fill {}{} buckets). Target bucket indices: {}..{}\n",
                LEVEL, ws.nb_buckets[LEVEL], LEVEL-1,
                LEVEL - 1, my_shorthint_t::rtti[0],
                first_region0_index,
                first_region0_index + ws.nb_buckets[LEVEL] * BRS[LEVEL] / BRS[1]);

        for (int i = 0; i < ws.nb_buckets[LEVEL]; ++i) {
            /* This is quite suspicious. Shouldn't we do BRS[LEVEL] /
             * BRS[LEVEL - 1] instead?
             */
            uint32_t const N = first_region0_index + i * (BRS[LEVEL] / BRS[1]);
            downsort_tree<LEVEL - 1>(ws, wc_p, aux_p, Q, pool, i, N,
                                     precomp_plattices, w);
        }
    } else {
        /* Now fill_in_buckets has completed for all levels. Time to check
         * that we had no overflow, and move on to process_bucket_region.
         */

        ws.check_buckets_max_full();

        if (print_slice_statistics)
            /* It may make sense to print a diagnosis of the fill ratio
             * of the different buckets. It's quite verbose, though */
            ws.slice_statistics(LEVEL);

        auto exc = pool.get_exceptions<buckets_are_full>(thread_pool::QUEUE_GENERIC);
        if (!exc.empty()) {
            throw *std::ranges::max_element(exc);
        }

        // it seems difficult to compute the max target bucket index, in
        // fact. Well of course it should be ws.nb_buckets[1], but just
        // based on the input that we have, it's less obvious.
        // size_t const(&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
        verbose_fmt_print(0, 3,
                "# calling process_bucket_region"
                " on regions of indices {}..\n",
                first_region0_index);

        /* PROCESS THE REGIONS AT LEVEL 0 */
        process_many_bucket_regions(ws, wc_p, aux_p, Q, pool, first_region0_index, w);

        /* We need that, because the next downsort_tree call in the loop
         * above (for LEVEL>1) will reset the pointers while filling the 1l
         * buckets -- and we read the 1l buckets from PBR.
         */
        if (ws.toplevel > 1)
            pool.drain_queue(thread_pool::QUEUE_GENERIC);
    }
}

template <int LEVEL>
void downsort_tree(
    nfs_work & ws,
    std::shared_ptr<nfs_work_cofac> wc_p,
    std::shared_ptr<nfs_aux> aux_p,
    ALGO::special_q_data const & Q,
    thread_pool & pool,
    uint32_t bucket_index, /* for the current level ! */
    uint32_t first_region0_index,
    std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> & precomp_plattices,
    where_am_I & w)
{
    if (ws.conf.needs_resieving()) {
        downsort_tree_inner<LEVEL, true>(ws, wc_p, aux_p, Q, pool, bucket_index,
                                         first_region0_index, precomp_plattices,
                                         w);
    } else {
        downsort_tree_inner<LEVEL, false>(ws, wc_p, aux_p, Q, pool, bucket_index,
                                          first_region0_index, precomp_plattices,
                                          w);
    }
}
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
