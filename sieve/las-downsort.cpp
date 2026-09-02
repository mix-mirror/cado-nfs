#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "fmt/format.h"

#include "bucket.hpp"
#include "las-auxiliary-data.hpp"
#include "las-bkmult.hpp"
#include "las-config.hpp"
#include "las-downsort.hpp"
#include "las-fill-in-buckets.hpp"


#include "las-threads.hpp"

#include "las-globals.hpp"
#include "las-process-bucket-region.hpp"
#include "las-report-stats.hpp"
#include "las-siever-config.hpp"
#include "las-smallsieve.hpp"   // IWYU pragma: keep
                                // SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE
#include "las-threads-work-data.hpp"
#include "las-where-am-i-proxy.hpp"
#include "macros.h"
#include "multityped_array.hpp"
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

template<bool WITH_HINTS>
struct downsort_object {
    nfs_work &ws;
    std::shared_ptr<nfs_work_cofac> wc_p;
    std::shared_ptr<nfs_aux> aux_p;
    ALGO::special_q_data const & Q;
    thread_pool &pool;
    std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> & precomp_plattices;
    where_am_I & w;
    using my_longhint_t = hints_proxy<WITH_HINTS>::l;
    using my_shorthint_t = hints_proxy<WITH_HINTS>::s;

    downsort_object(
        nfs_work &ws,
        std::shared_ptr<nfs_work_cofac> const & wc_p,
        std::shared_ptr<nfs_aux> const & aux_p,
        ALGO::special_q_data const & Q,
        thread_pool &pool,
        std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> & precomp_plattices,
        where_am_I & w)
        : ws(ws)
        , wc_p(wc_p)
        , aux_p(aux_p)
        , Q(Q)
        , pool(pool)
        , precomp_plattices(precomp_plattices)
        , w(w)
    {}

    /* This is auxiliary only. We downsort stuff that we already downsorted.
     * So it applies only if LEVEL+1 is itself not the toplevel.
     * For this reason, we must have a specific instantiation that reduces
     * this to a no-op if LEVEL+1>=3, because there's no longhint_t for level
     * 3 presently.
     */
    template <int LEVEL>
    void ds_aux(task_group & tg,
            int side,
            uint32_t bucket_index)
    {
        static_assert(LEVEL <= MAX_TOPLEVEL - 1);

        nfs_work::side_data & wss(ws.sides[side]);

        // Early exit if this side does not reach LEVEL + 2, since for
        // LEVEL+1 to contain longhints, we need the toplevel to be LEVEL+2
        // or above.
        if (wss.fbs->get_toplevel() < LEVEL + 2)
            return;

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
                [this, side, bucket_index, &BA_in](worker_thread * worker, where_am_I && w) {
                    nfs_work::side_data & wss(ws.sides[side]);
                    fb_factorbase::slicing const & fbs(*wss.fbs);
                    int const id = worker->rank();
                    nfs_aux::thread_data & taux(aux_p->th[id]);
                    taux.w = std::move(w);
                    timetree_t & timer(aux_p->get_timer(worker));
                    ENTER_THREAD_TIMER(timer);
                    MARK_TIMER_FOR_SIDE(timer, side);
                    taux.w = std::move(w);
                    CHILD_TIMER(timer, fmt::format("downsort<{}>", LEVEL));
                    auto tt = worker->trace(chronograms::DS(side,LEVEL,bucket_index));
                    downsort<LEVEL + 1>(fbs,
                            wss.acquire_BA<LEVEL, my_longhint_t>(wss.rank_BA(BA_in)),
                            BA_in, bucket_index, taux.w);
                },
                where_am_I(w));
        }
    }
    // }}}
    // {{{ FIB (fill-in-buckets)
    // For internal levels, the fill-in is not exactly the same as for
    // top-level, since the plattices have already been precomputed.
    template <int LEVEL>
    void fib_internal(
            worker_thread * worker,
            int side,
            plattices_vector_t & plattices_vector,
            uint32_t const first_region0_index)
    {
        /* Import some contextual stuff */
        int const id = worker->rank();
        nfs_aux::thread_data & taux(aux_p->th[id]);
        timetree_t & timer(aux_p->get_timer(worker));
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
        WHERE_AM_I_UPDATE(w, i, plattices_vector.get_index());
        WHERE_AM_I_UPDATE(w, N, first_region0_index);

        try {
            auto acquired = wss.reserve_BA<LEVEL, my_shorthint_t>();
            auto tt = worker->trace(chronograms::FIB(
                        side,
                        LEVEL,
                        wss.rank_BA(acquired.access()),
                        plattices_vector.get_index()));
                
            /* Get an unused bucket array that we can write to */
            /* clearly, reserve_BA() possibly throws. As it turns out,
             * fill_in_buckets_lowlevel<> does not, at least currently. One
             * could imagine that it could throw, so let's wrap it too.
             */
            fill_in_buckets_lowlevel<LEVEL, my_shorthint_t>(
                    acquired.access(),
                    ws, Q, plattices_vector,
                    first_region0_index, w);
        } catch (buckets_are_full & e) {
            e.side = side;
            throw e;
        }
    }

    // {{{ DS (downsort) -- should refactor shorthint and longhint.
    template <int LEVEL>
        void ds(task_group & tg, int side, uint32_t bucket_index)
    {
        /* FIRST: Downsort what is coming from the level above, for this
         * bucket index */
        // All these BA are global stuff; see reservation_group.
        // We reserve those where we write, and access the ones for
        // reading without reserving. We require that things at level
        // above are finished before entering here.

        nfs_work::side_data & wss(ws.sides[side]);

        if (wss.fbs->get_toplevel() < LEVEL + 1)
            return;

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
                    tg,
                    [this, side, &BA_in, bucket_index](worker_thread * worker, where_am_I && w) {
                    nfs_work::side_data & wss(ws.sides[side]);
                    fb_factorbase::slicing const & fbs(*wss.fbs);
                    int const id = worker->rank();
                    nfs_aux::thread_data & taux(aux_p->th[id]);
                    timetree_t & timer(aux_p->get_timer(worker));
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


    // first_region0_index is a way to remember where we are in the tree.
    // The depth-first is a way to process all the the regions of level 0 in
    // increasing order of j-value.
    // first_region0_index * nb_lines_per_region0 therefore gives the j-line
    // where we are. This is what is called N by WHERE_AM_I and friends.
    template<int LEVEL>
    void fib(task_group & tg,
            int side,
            uint32_t first_region0_index)
    {
        /* SECOND: fill in buckets at this level, for this region. */
        nfs_work::side_data & wss(ws.sides[side]);
        wss.reset_all_pointers<LEVEL, my_shorthint_t>();

        auto & BA_outs = wss.bucket_arrays<LEVEL, my_shorthint_t>();
        auto & lattices = precomp_plattices[side].template get<LEVEL>();

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
            pool.add_task(tg, thread_pool::QUEUE_GENERIC, it.get_weight(),
                    [this, side, &it, first_region0_index](worker_thread * worker, where_am_I && w) {
                    int const id = worker->rank();
                    nfs_aux::thread_data & taux(aux_p->th[id]);
                    taux.w = std::move(w);
                    fib_internal<LEVEL>(worker, side, it, first_region0_index);
                    },
                    where_am_I(w));
        }
    }

    // {{{ SSS (start small sieve), only for level 1.
    void sss(task_group & tg, int side, uint32_t first_region0_index)
    {
        /* Prepare for PBR: we need to precompute the small sieve positions
         * for all the small sieved primes.
         *
         * For ws.toplevel==1, we don't reach here, of course, and the
         * corresponding initialization is done with identical code in
         * las.cpp
         */
        ASSERT(ws.toplevel > 1);
        nfs_work::side_data const & wss(ws.sides[side]);

        wss.ssd->small_sieve_prepare_many_start_positions(
                pool, &tg,
                first_region0_index,
                std::min(SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE, ws.nb_buckets[1]),
                ws.conf.logI, Q.sublat);
        tg.on_complete([&wss]() {
                wss.ssd->small_sieve_activate_many_start_positions();
                });
    }


    template <int LEVEL>
    void fib_ds_sss(
            uint32_t bucket_index, /* for the current level ! */
            uint32_t first_region0_index
            )
    {
        /* LEVEL is not the toplevel here, so we must have the following: */
        static_assert(LEVEL <= MAX_TOPLEVEL - 1);

        int const nsides = ws.las.cpoly.nsides();
        nfs_aux & aux(*aux_p);
        timetree_t & timer(aux.rt.timer);


        CHILD_TIMER(timer, fmt::format("downsort_tree<{}>", LEVEL));
        TIMER_CATEGORY(timer, sieving_mixed());
        ASSERT_ALWAYS(LEVEL > 0);

        WHERE_AM_I_UPDATE(w, N, first_region0_index);

        std::vector<task_group> ds_tgs(nsides);
        std::vector<task_group> ds_aux_tgs(nsides);
        std::vector<task_group> fib_tgs(nsides);
        std::vector<task_group> sss_tgs(nsides);

        for (int side = 0; side < nsides; ++side) {
            nfs_work::side_data const & wss(ws.sides[side]);
            if (wss.no_fb())
                continue;

            WHERE_AM_I_UPDATE(w, side, side);
            TIMER_CATEGORY(timer, sieving(side));

            auto & ds_tg = ds_tgs[side];
            auto & ds_aux_tg = ds_aux_tgs[side];
            auto & fib_tg = fib_tgs[side];
            auto & sss_tg = sss_tgs[side];

            ds<LEVEL>(ds_tg, side, bucket_index);

            if (LEVEL < ws.toplevel - 1) {
                /* ds_aux() (the "downsort what was already downsorted"
                 * pass) appends to the very same <LEVEL, my_longhint_t>
                 * bucket arrays that ds() fills from the shorthint
                 * updates: ds() creates the single slice with
                 * add_slice_index(0), ds_aux() only appends. So ds_aux()
                 * must run only once ds_tg (the shorthint downsort) has
                 * completed for this side -- never concurrently. */
                ds_tg.on_complete(
                        [this, &ds_aux_tg, side, bucket_index]() mutable {
                            ds_aux<LEVEL>(ds_aux_tg, side, bucket_index);
                        });
            }

            fib<LEVEL>(fib_tg, side, first_region0_index);

            if (LEVEL == 1) {
                ASSERT(ws.toplevel > 1);

                sss(sss_tg, side, first_region0_index);
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
    }

    template<int LEVEL> void recurse(
            uint32_t first_region0_index)
    {
        static_assert(LEVEL>1);
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
            tree<LEVEL - 1>(i, N);
        }
    }

    void pbr(uint32_t first_region0_index)
    {
        /* Now fill_in_buckets has completed for all levels. Time to check
         * that we had no overflow, and move on to process_bucket_region.
         */

        ws.check_buckets_max_full();

        if (print_slice_statistics)
            /* It may make sense to print a diagnosis of the fill ratio
             * of the different buckets. It's quite verbose, though */
            /* XXX should this go in tree<> instead? */
            ws.slice_statistics(1);

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

    template<int LEVEL> void tree(
            uint32_t bucket_index,
            uint32_t first_region0_index)
    {
        if constexpr (LEVEL > 0)
            fib_ds_sss<LEVEL>(bucket_index, first_region0_index);
        if constexpr (LEVEL > 1) {
            /* Because of this "if constexpr", the recursion terminates */
            recurse<LEVEL>(first_region0_index);
        } else {
            pbr(first_region0_index);
        }
    }

    void tree_toplevel()
    {
        // Visit the downsorting tree depth-first.
        // If toplevel = 1, then this is just processing all bucket
        // regions.
        size_t const(&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;

        static_assert(MAX_TOPLEVEL <= 3);

        /* TODO: it's ugly. */

        for (int i = 0; i < ws.nb_buckets[ws.toplevel]; i++) {
            if (ws.task->must_take_decision())
                break;
            switch (ws.toplevel) {
                case 1:
                    /* there should only be one loop, then. */
                    ASSERT_ALWAYS(i == 0);
                    tree<0>(0, 0);
                    return;
#if MAX_TOPLEVEL >= 2
                case 2:
                    tree<1>(i, i*BRS[2]/BRS[1]);
                    break;
#endif
#if MAX_TOPLEVEL >= 3
                case 3:
                    tree<2>(i, i*BRS[3]/BRS[1]);
                    break;
#endif
                default:
                    ASSERT_ALWAYS(0);
            }
        }
    }
};

void downsort_toplevel(
        nfs_work &ws,
        std::shared_ptr<nfs_work_cofac> const & wc_p,
        std::shared_ptr<nfs_aux> const & aux_p,
        ALGO::special_q_data const & Q,
        thread_pool &pool,
        std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> & precomp_plattices,
        where_am_I & w)
{
    if (ws.conf.needs_resieving()) {
        downsort_object<true> D(ws, wc_p, aux_p, Q, pool, precomp_plattices, w);
        D.tree_toplevel();
    } else {
        downsort_object<false> D(ws, wc_p, aux_p, Q, pool, precomp_plattices, w);
        D.tree_toplevel();
    }
}

// some explicit instantiations are needed to terminate the compile-time
// recursions. The code is such that at runtime, we never reach here!
template <>
template <>
void downsort_object<true>::ds_aux<2>(task_group &, int, uint32_t) { ASSERT_ALWAYS(0); }
template <>
template <>
void downsort_object<false>::ds_aux<2>(task_group &, int, uint32_t) { ASSERT_ALWAYS(0); }
