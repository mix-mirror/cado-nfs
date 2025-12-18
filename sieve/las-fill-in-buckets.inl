#ifndef CADO_SIQS_FILL_IN_BUCKETS_INL
#define CADO_SIQS_FILL_IN_BUCKETS_INL

/***********************************************************************/
/* multithreaded processing of make_lattice_bases (a.k.a
 * precomp_plattices)
 *
 * This creates one control object per slice, with storage ownership of
 * the control object transfered to the called function. Because we
 * depend on the slice, the type of the object is parameterized by the
 * slice type.
 *
 * We may elect to make the "model" a shared_ptr.
 */

template <int LEVEL>
struct make_lattice_bases_parameters_base : public task_parameters {
    int side;
    nfs_work & ws;
    qlattice_basis const & Q;
    precomp_plattice_t<LEVEL> & V;
    make_lattice_bases_parameters_base(int side, nfs_work & ws,
            qlattice_basis const & Q,
            precomp_plattice_t<LEVEL> & V)
        : side(side)
        , ws(ws)
        , Q(Q)
        , V(V)
    {
    }
};
template <int LEVEL, class FB_ENTRY_TYPE>
struct make_lattice_bases_parameters
    : public make_lattice_bases_parameters_base<LEVEL> {
    using super = make_lattice_bases_parameters_base<LEVEL>;
    fb_slice<FB_ENTRY_TYPE> const & slice;
    make_lattice_bases_parameters(super const & model,
                                  fb_slice<FB_ENTRY_TYPE> const & slice)
        : super(model)
        , slice(slice)
    {
    }
};

template <int LEVEL, class FB_ENTRY_TYPE>
static task_result * make_lattice_bases(worker_thread * worker MAYBE_UNUSED,
                                        task_parameters * _param, int)
{
    auto const * param =
        static_cast<
            make_lattice_bases_parameters<LEVEL, FB_ENTRY_TYPE> const *>(
            _param);

    nfs_work const & ws(param->ws);
    auto const & Q(param->Q);
    int const logI = ws.conf.logI;
    sublat_t const & sublat(Q.sublat);
    auto & V(param->V);
    auto const & slice(param->slice);

    auto const index0 = ws.sides[param->side].fbs->get_part(LEVEL).first_slice_index;
    auto const index = slice.get_index();
    auto const relative_index = index - index0;
    ASSERT_ALWAYS(relative_index < V.size());

    typename FB_ENTRY_TYPE::transformed_entry_t transformed;
    /* Create a transformed vector and store the index of the fb_slice we
     * currently transform */

    /* we don't really need the fence at this point, except that we do one
     * shot of "next()", and that happens to need the fence. Only logI is
     * needed, though.
     */
    plattice_enumerator::fence const F(ws.conf.logI, 0);

    plattices_vector_t result(index, slice.get_weight());
    slice_offset_t i_entry = 0;
    for (auto const & e: slice) {
        increment_counter_on_dtor<slice_offset_t> const _dummy(i_entry);
        if (!Q.is_coprime_to(e.p))
            continue;
        if (discard_power_for_bucket_sieving(e))
            continue;
        e.transform_roots(transformed, Q);
        for (unsigned char i_root = 0; i_root != transformed.nr_roots;
             i_root++) {
            fbroot_t const r = transformed.get_r(i_root);
            bool const proj = transformed.get_proj(i_root);
            plattice_info const pli =
                plattice_info(transformed.get_q(), r, proj, logI);
            plattice_enumerator ple(pli, i_entry, logI, sublat);
            // Skip (0,0) unless we have sublattices.
            if (!sublat.m)
                ple.next(F);
            if (LIKELY(!pli.is_discarded()))
                result.push_back(ple);
        }
    }
    /* This is moved, not copied. Note that V is a reference. */
    V[relative_index] = std::move(result);
    delete param;
    return new task_result;
}

void fill_in_buckets_prepare_plattices(
    nfs_work & ws, qlattice_basis const & Q, thread_pool & pool, int side,
    cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS-1> & precomp_plattice)
{
    /* this will *not* do anything for level==ws.toplevel, by design */
    precomp_plattice.foreach([&](auto & precomp_plattice) {
        /* T is precomp_plattice_t<n> for some level n */
        using T = std::remove_reference_t<decltype(precomp_plattice)>;
        if (T::level >= ws.toplevel)
            return;

        nfs_work::side_data const & wss(ws.sides[side]);
        fb_factorbase::slicing::part const & P = wss.fbs->get_part(T::level);
        /* We pre-assign the result, so that all threads can write to it
         * comfortably.
         *
         * It would be nice to have a way to notify that all threads here are
         * done with their job.
         */
        precomp_plattice.assign(P.nslices(), plattices_vector_t());
        make_lattice_bases_parameters_base<T::level> const model {side, ws, Q, precomp_plattice};
        P.slices.foreach([&](auto const & sl) {
                for(auto const & s : sl) {
                    using E = std::remove_reference_t<decltype(s)>::entry_t;
                    auto param = new make_lattice_bases_parameters<T::level, E>(model, s);
                    task_function_t f = make_lattice_bases<T::level, E>;
                    pool.add_task(f, param, 0);
                }
        });
    });
}

/***********************************************************************/
/* multithreaded processing of fill_in_buckets_toplevel (both with and
 * without sublattices) is more complicated. First because the important
 * functions are not the ones whose prototype is the one we expect most
 * from a multithreade task, second because we strive to manage
 * exceptions properly. So we go through several quirky paths below.
 */

// At top level, the fill-in of the buckets must interleave
// the root transforms and the FK walks, otherwise we spend a lot of time
// doing nothing while waiting for memory.
//
// With Sublat, this function can have two modes:
//   - process the given slice, and store the corresponding FK-basis in
//     precomp_slice for later use.
//   - use the pre-processed precomputed FK_basis.
// If (and only if) we are dealing with (i,j) == (0,1) mod m,
// we are in the second mode.
//
//
// FIXME FIXME FIXME: tons of duplicated code, here!!!
// But putting if() in critical loops can kill performance (I tried...)

template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static void fill_in_buckets_toplevel_sublat(
    bucket_array_t<LEVEL, TARGET_HINT> & orig_BA, nfs_work & ws,
    qlattice_basis const & Q,
    fb_slice<FB_ENTRY_TYPE> const & slice,
    plattices_dense_vector_t * p_precomp_slice, where_am_I & w)
{
    int const logI = ws.conf.logI;

    plattices_dense_vector_t & precomp_slice(*p_precomp_slice);

    ASSERT_ALWAYS(Q.sublat.m);
    bool const first_sublat = Q.sublat.i0 == 0 && Q.sublat.j0 == 1;
    /* local copy. Gain a register + use stack */
    bucket_array_t<LEVEL, TARGET_HINT> BA = std::move(orig_BA);

    slice_index_t const slice_index = slice.get_index();

    /* Write new set of pointers for the new slice */
    BA.add_slice_index(slice_index);

    typename FB_ENTRY_TYPE::transformed_entry_t transformed;

    /* top level: the fence we care about is the one defined by J */
    plattice_enumerator::fence const F(ws.conf.logI, ws.J);

    int logB = LOG_BUCKET_REGIONS[LEVEL];
    typename bucket_array_t<LEVEL, TARGET_HINT>::update_t::br_index_t bmask =
        (1UL << logB) - 1;

    // FIXME: A LOT OF DUPLICATED CODE, HERE!!!
    if (first_sublat) {
        slice_offset_t i_entry = 0;
        for (auto const & e: slice) {
            increment_counter_on_dtor<slice_offset_t> const _dummy(i_entry);
            if (!Q.is_coprime_to(e.p))
                continue;
#ifdef BUCKET_SIEVE_POWERS
            /* the combination of bucket-sieving powers + sublattices means
             * that powers of the primes that divide the sublattice determinant
             * may be bucket-sieved. And of course, that leads to problems.
             *
             * technically, Q.sublat.m could be composite, in which case we
             * would have a gcd to compute, here. The only really useful case
             * at the moment is m=3 though.
             */
            if (Q.sublat.m == 3) {
                if (e.p == 3)
                    continue;
            } else if (gcd_ul(e.p, Q.sublat.m) > 1) {
                continue;
            }
#endif
            if (discard_power_for_bucket_sieving(e))
                continue;
            e.transform_roots(transformed, Q);
            for (unsigned char i_root = 0; i_root != transformed.nr_roots;
                 i_root++) {
                fbroot_t const r = transformed.get_r(i_root);
                bool const proj = transformed.get_proj(i_root);
                plattice_info const pli(transformed.get_q(), r, proj, logI);
                // In sublat mode, save it for later use
                precomp_slice.push_back(plattice_info_dense_t(pli, i_entry));

                plattice_enumerator ple(pli, i_entry, logI, Q.sublat);

                if (ple.done(F))
                    continue;
                if (pli.is_discarded())
                    continue;
                slice_offset_t const hint = ple.get_hint();
                ASSERT(hint == i_entry);
                WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
                const fbprime_t p = slice.get_prime(hint);
                WHERE_AM_I_UPDATE(w, p, p);
#else
                const fbprime_t p = 0;
#endif

                typename bucket_array_t<LEVEL, TARGET_HINT>::update_t u(
                    0, p, hint, slice_index);

                // Handle the rare special cases
                // XXX Here, we're not bucket-sieving projective primes at
                // all, and neither do we bucket-sieve primes with root equal
                // to zero.
                if (UNLIKELY(pli.is_vertical_line(logI) ||
                             pli.is_projective_like(logI)))
                    continue;

                /* Now, do the real work: the filling of the buckets */
                // Without sublattices, we test (very basic) coprimality,
                while (!ple.done(F)) {
                    u.set_x(ple.get_x() & bmask);
                    BA.push_update(ple.get_x() >> logB, u, w);
                    ple.next(F);
                }
            }
        }
    } else { // Use precomputed FK-basis
        for (unsigned int i = 0; i < precomp_slice.size(); ++i) {
            plattice_info const pli(precomp_slice[i].unpack(logI));
            slice_offset_t const i_entry = precomp_slice[i].get_hint();

            plattice_enumerator ple(pli, i_entry, logI, Q.sublat);

            if (ple.done(F))
                continue;
            if (pli.is_discarded())
                continue;
            slice_offset_t const hint = ple.get_hint();
            WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
            const fbprime_t p = slice.get_prime(hint);
            WHERE_AM_I_UPDATE(w, p, p);
#else
            const fbprime_t p = 0;
#endif

            typename bucket_array_t<LEVEL, TARGET_HINT>::update_t u(
                0, p, hint, slice_index);

            // Handle (well, do not handle, in fact) the rare special cases
            if (UNLIKELY(pli.is_vertical_line(logI) ||
                         pli.is_projective_like(logI)))
                continue;

            /* Now, do the real work: the filling of the buckets */
            // Without sublattices, we test (very basic) coprimality,
            // otherwise not atm. FIXME!
            while (!ple.done(F)) {
                u.set_x(ple.get_x() & bmask);
                BA.push_update(ple.get_x() >> logB, u, w);
                ple.next(F);
            }
        }
    }
    // printf("%.3f\n", BA.max_full());
    orig_BA = std::move(BA);
}

/* TARGET_HINT is shorthint_t or void */
template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static void
fill_in_buckets_toplevel(bucket_array_t<LEVEL, TARGET_HINT> & orig_BA,
                         nfs_work & ws, fb_slice<FB_ENTRY_TYPE> const & slice,
                         qlattice_basis const & Q,
                         plattices_dense_vector_t * /* unused */,
                         where_am_I & w)
{
    int const logI = ws.conf.logI;

    ASSERT_ALWAYS(!Q.sublat.m);

    /* local copy. Gain a register + use stack */
    bucket_array_t<LEVEL, TARGET_HINT> BA = std::move(orig_BA);

    slice_index_t const slice_index = slice.get_index();

    /* Write new set of pointers for the new slice */
    BA.add_slice_index(slice_index);
    WHERE_AM_I_UPDATE(w, i, slice_index);

    typename FB_ENTRY_TYPE::transformed_entry_t transformed;

    /* top level: the fence we care about is the one defined by J */
    plattice_enumerator::fence const F(ws.conf.logI, ws.J);

    slice_offset_t i_entry = 0;

    /* yes, we want the level-1 regions here */
    int logB1 = LOG_BUCKET_REGIONS[1];
    uint32_t maskB1I = (UINT32_C(1) << std::min(logB1, logI)) - 1;

    int logB = LOG_BUCKET_REGIONS[LEVEL];
    typename bucket_array_t<LEVEL, TARGET_HINT>::update_t::br_index_t bmask =
        (1UL << logB) - 1;

    for (auto const & e: slice) {
        increment_counter_on_dtor<slice_offset_t> const _dummy(i_entry);
        if (!Q.is_coprime_to(e.p))
            continue;
        if (discard_power_for_bucket_sieving(e))
            continue;
        e.transform_roots(transformed, Q);
        for (unsigned char i_root = 0; i_root != transformed.nr_roots;
             i_root++) {
            fbroot_t const r = transformed.get_r(i_root);
            bool const proj = transformed.get_proj(i_root);
            plattice_info const pli(transformed.get_q(), r, proj, logI);

            plattice_enumerator ple(pli, i_entry, logI);

            // Skip (i,j)=(0,0)
            ple.next(F);

            if (pli.is_discarded())
                continue;

            slice_offset_t const hint = ple.get_hint();
            WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
            const fbprime_t p = slice.get_prime(hint);
            WHERE_AM_I_UPDATE(w, p, p);
#else
            const fbprime_t p = 0;
#endif

            typename bucket_array_t<LEVEL, TARGET_HINT>::update_t u(
                0, p, hint, slice_index);

            // Handle the rare special cases
            /* projective-like:
             *
             * ple sets its first position in the (i,j) plane to (1,0),
             * which will typically be the _only_ hit in the normal case.
             *
             * there are more subtle cases that can show up though, because
             * of projective powers, and the combination with
             * adjust-strategy 2 (see bug 30012).
             *
             * the first hit (and only hit on the first line) can be (g,0)
             * for any g. but other lines may hit.
             */
            /* vertical:
             *
             * Root=0: only update is at (0,something).
             * note that "something" might be large !
             */
            if (UNLIKELY(ple.is_projective_like(logI))) {
                while (!ple.done(F)) {
                    u.set_x(ple.get_x() & bmask);
                    int N = ple.get_x() >> logB;
                    int n =
                        ple.advance_to_end_of_row_or_smallest_region(maskB1I);
                    BA.push_row_update(slice_index, ple.get_inc_step(),
                            N, n, u, w);
                    ple.next(F);
                }
            } else if (UNLIKELY(pli.is_vertical_line(logI))) {
                if (!ple.done(F)) {
                    u.set_x(ple.get_x() & bmask);
                    BA.push_update(ple.get_x() >> logB, u, w);
                    ple.finish();
                }
            } else {
                /* Now, do the real work: the filling of the buckets */
                while (!ple.done(F)) {
                    if (LIKELY(ple.probably_coprime(F))) {
                        u.set_x(ple.get_x() & bmask);
                        BA.push_update(ple.get_x() >> logB, u, w);
                    }
                    ple.next(F);
                }
            }
        }
    }
    // printf("%.3f\n", BA.max_full());
    orig_BA = std::move(BA);
}

/* TARGET_HINT is shorthint_t or void */
template <int LEVEL, typename TARGET_HINT>
static void
fill_in_buckets_lowlevel(bucket_array_t<LEVEL, TARGET_HINT> & orig_BA,
                         nfs_work & ws,
                         qlattice_basis const & Q,
                         plattices_vector_t & plattices_vector,
                         bool first_reg MAYBE_UNUSED, where_am_I & w)
{
    fmt::println(stderr, "#DEV {} LEVEL={} J={} logB={} logI={}",
                         __func__, LEVEL, ws.J, LOG_BUCKET_REGIONS[LEVEL],
                         ws.conf.logI);
    int const logI = ws.conf.logI;

    /* The timer stuff is dealt with by the caller */
    slice_index_t const slice_index = plattices_vector.get_index();

    /* local copy. Gain a register + use stack */
    bucket_array_t<LEVEL, TARGET_HINT> BA = std::move(orig_BA);

    /* Write new set of pointers for the new slice */
    BA.add_slice_index(slice_index);
    WHERE_AM_I_UPDATE(w, i, slice_index);

    /* we used to look up BUCKET_REGIONS[LEVEL + 1] here, which doesn't
     * really seem to make sense. I expect that ws.J actually yields a
     * stricter bound in all cases.
     */
    plattice_enumerator::fence const F(ws.conf.logI, ws.J,
            (LEVEL + 1 < FB_MAX_PARTS ? BUCKET_REGIONS[LEVEL + 1] : SIZE_MAX));
    /* just checking... */
    if (LEVEL == FB_MAX_PARTS - 1)
        ASSERT_ALWAYS((ws.J << ws.conf.logI) < (BUCKET_REGIONS[LEVEL] << 8));

    /* yes, we want the level-1 regions here */
    int logB1 = LOG_BUCKET_REGIONS[1];
    uint32_t maskB1I = (UINT32_C(1) << std::min(logB1, logI)) - 1;

    int logB = LOG_BUCKET_REGIONS[LEVEL];
    typename bucket_array_t<LEVEL, TARGET_HINT>::update_t::br_index_t bmask =
        (1UL << logB) - 1;

    for (auto & ple_orig: plattices_vector) {
        // Work with a copy, otherwise we don't get all optimizations.
        // Maybe with a wise use of the 'restrict' keyword, we might get
        // what we want, but this is C++11, anyway.
        //
        // FIXME we're c++11 now. Look into this.
        plattice_enumerator ple(ple_orig);

        slice_offset_t const hint = ple.get_hint();
        WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
        /* this is a bit expensive, since we're scanning all parts.
         * Fortunately it's only a debug call anyway. */
        fb_slice_interface const & slice =
            (*w->sides[w->side].fbs)[slice_index];
        fbprime_t const p = slice.get_prime(hint);
        WHERE_AM_I_UPDATE(w, p, p);
#else
        const fbprime_t p = 0;
#endif

        typename bucket_array_t<LEVEL, TARGET_HINT>::update_t u(0, p, hint,
                                                                slice_index);

        // Handle the rare special cases
        /* see fill_in_bucket_toplevel. */
        if (UNLIKELY(ple.is_projective_like(logI))) {
            if (Q.sublat.m)
                continue; /* XXX headaches ! */

            while (!ple.done(F)) {
                u.set_x(ple.get_x() & bmask);
                int N = ple.get_x() >> logB;
                int n = ple.advance_to_end_of_row_or_smallest_region(maskB1I);
                BA.push_row_update(slice_index, ple.get_inc_step(), N, n, u, w);
                ple.next(F);
            }
            /* we now do the end of loop normally: store x into ple_orig, and
             * then advance to the next area. This is because more rows can
             * be interesting as we go towards increasing j's
             */
        } else if (UNLIKELY(ple.is_vertical_line(logI))) {
            if (Q.sublat.m)
                continue; /* XXX headaches ! */

            if (!ple.done(F)) {
                u.set_x(ple.get_x() & bmask);
                BA.push_update(ple.get_x() >> logB, u, w);
                // ple.next(F);
                ple.finish();
            }
        } else {
            /* Now, do the real work: the filling of the buckets */
            // Without sublattices, we test (very basic) coprimality,
            // otherwise not atm. FIXME!
            if (!Q.sublat.m) {
                while (!ple.done(F)) {
                    if (LIKELY(ple.probably_coprime(F))) {
                        u.set_x(ple.get_x() & bmask);
                        BA.push_update(ple.get_x() >> logB, u, w);
                    }
                    ple.next(F);
                }
            } else {
                while (!ple.done(F)) {
                    u.set_x(ple.get_x() & bmask);
                    BA.push_update(ple.get_x() >> logB, u, w);
                    ple.next(F);
                }
            }
        }
        // save current position, and prepare for next area.
        ple_orig.set_x(ple.get_x());
        ple_orig.advance_to_next_area(F);
    }
    // printf("%.3f\n", BA.max_full());
    orig_BA = std::move(BA);
}

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
static void downsort_aux(fb_factorbase::slicing const & fbs, nfs_work & ws,
                         nfs_aux & aux, thread_pool & pool, int side,
                         uint32_t bucket_index, where_am_I & w)
{
    static_assert(LEVEL <= MAX_TOPLEVEL - 1);

    using my_longhint_t = hints_proxy<WITH_HINTS>::l;

    nfs_work::side_data & wss(ws.sides[side]);

    auto const & BA_ins = wss.bucket_arrays<LEVEL + 1, my_longhint_t>();
    auto & BA_outs = wss.bucket_arrays<LEVEL, my_longhint_t>();

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
        pool.add_task_lambda(
            [&, side, w](worker_thread * worker, int bucket_index) {
                nfs_aux::thread_data & taux(aux.th[worker->rank()]);
                timetree_t & timer(aux.get_timer(worker));
                ENTER_THREAD_TIMER(timer);
                MARK_TIMER_FOR_SIDE(timer, side);
                taux.w = w;
                CHILD_TIMER(timer, TEMPLATE_INST_NAME(downsort, LEVEL));
                auto & BA_out(
                    wss.reserve_BA<LEVEL, my_longhint_t>(wss.rank_BA(BA_in)));
                downsort<LEVEL + 1>(fbs, BA_out, BA_in, bucket_index, taux.w);
                wss.template release_BA<LEVEL, my_longhint_t>(BA_out);
            },
            bucket_index, 0);
    }
}

#if MAX_TOPLEVEL == 2
template <>
void downsort_aux<1, false>(fb_factorbase::slicing const &, nfs_work &, nfs_aux &,
                     thread_pool &, int, uint32_t, where_am_I &)
{
}
template <>
void downsort_aux<1, true>(fb_factorbase::slicing const &, nfs_work &, nfs_aux &,
                     thread_pool &, int, uint32_t, where_am_I &)
{
}
#endif
#if MAX_TOPLEVEL == 3
template <>
void downsort_aux<2, false>(fb_factorbase::slicing const &, nfs_work &, nfs_aux &,
                     thread_pool &, int, uint32_t, where_am_I &)
{
}
template <>
void downsort_aux<2, true>(fb_factorbase::slicing const &, nfs_work &, nfs_aux &,
                     thread_pool &, int, uint32_t, where_am_I &)
{
}
#endif
static_assert(MAX_TOPLEVEL == 3);

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
    qlattice_basis const & Q,
    thread_pool & pool,
    uint32_t bucket_index, /* for the current level ! */
    uint32_t first_region0_index,
    std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> & precomp_plattices,
    where_am_I & w)
{
    /* LEVEL is not the toplevel here, so we must have the following: */
    static_assert(LEVEL <= MAX_TOPLEVEL - 1);

    int const nsides = ws.las.cpoly->nb_polys;
    nfs_aux & aux(*aux_p);
    timetree_t & timer(aux.rt.timer);

    using my_longhint_t = hints_proxy<WITH_HINTS>::l;
    using my_shorthint_t = hints_proxy<WITH_HINTS>::s;

    CHILD_TIMER(timer, TEMPLATE_INST_NAME(downsort_tree, LEVEL));
    TIMER_CATEGORY(timer, sieving_mixed());
    ASSERT_ALWAYS(LEVEL > 0);

    WHERE_AM_I_UPDATE(w, N, first_region0_index);

    for (int side = 0; side < nsides; ++side) {
        nfs_work::side_data & wss(ws.sides[side]);
        if (wss.no_fb())
            continue;

        WHERE_AM_I_UPDATE(w, side, side);
        TIMER_CATEGORY(timer, sieving(side));
        /* FIRST: Downsort what is coming from the level above, for this
         * bucket index */
        // All these BA are global stuff; see reservation_group.
        // We reserve those where we write, and access the ones for
        // reading without reserving. We require that things at level
        // above are finished before entering here.

        {
            /* This is the "dictionary" that maps slice indices to actual fb
             * entries. We rarely need it, except when downsorting short entries
             * in the case where we've eliminated the hint
             */
            fb_factorbase::slicing const & fbs(*wss.fbs);

            auto const & BA_ins = wss.bucket_arrays<LEVEL + 1, my_shorthint_t>();
            auto & BA_outs = wss.bucket_arrays<LEVEL, my_longhint_t>();

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

            /* otherwise the code here can't work */
            ASSERT_ALWAYS(BA_ins.size() == BA_outs.size());

            /* We create one output array for each input array, and we
             * process them in parallel. There would be various ways to
             * achieve that.
             */
            for (auto const & BA_in: BA_ins) {
                pool.add_task_lambda(
                    [&, side, w](worker_thread * worker, int bucket_index) {
                        nfs_aux::thread_data & taux(aux.th[worker->rank()]);
                        timetree_t & timer(aux.get_timer(worker));
                        taux.w = w;
                        ENTER_THREAD_TIMER(timer);
                        MARK_TIMER_FOR_SIDE(timer, side);
                        CHILD_TIMER(timer, TEMPLATE_INST_NAME(downsort, LEVEL));
                        auto & BA_out(wss.reserve_BA<LEVEL, my_longhint_t>(
                            wss.rank_BA(BA_in)));
                        BA_out.reset_pointers();
                        downsort<LEVEL + 1>(fbs, BA_out, BA_in, bucket_index,
                                            taux.w);
                        //wss.template release_BA<LEVEL, my_longhint_t>(BA_out);
                        wss.release_BA(BA_out);
                    },
                    bucket_index, 0);
            }
            // What comes from already downsorted data above. We put this in
            // an external function because we need the code to be elided or
            // LEVEL >= 2.
            if (LEVEL < ws.toplevel - 1) {
                pool.drain_queue(0);
                downsort_aux<LEVEL, WITH_HINTS>(fbs, ws, aux, pool, side, bucket_index, w);
            }
        }

        /* There might be a performance hit here, and honestly I'm not
         * 100% sure it's useful. The F9_sieve_3_levels test wants it,
         * and apparently really wants it here.
         */
        pool.drain_queue(0);

        {
            /* SECOND: fill in buckets at this level, for this region. */
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
                pool.add_task(
                        fill_in_buckets_one_slice_internal<LEVEL, my_shorthint_t>,
                        new fill_in_buckets_parameters<LEVEL> {
                        ws, aux, Q, side, (fb_slice_interface *)NULL, &it, NULL,
                        first_region0_index, const_ref(w)},
                        0, 0, it.get_weight());
            }
        }
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
            pool.add_task_lambda(
                [=, &ws, &aux](worker_thread * worker, int) {
                    timetree_t & timer(aux.get_timer(worker));
                    ENTER_THREAD_TIMER(timer);
                    MARK_TIMER_FOR_SIDE(timer, side);
                    SIBLING_TIMER(timer, "prepare small sieve");
                    nfs_work::side_data & wss(ws.sides[side]);
                    // if (wss.no_fb()) return;
                    SIBLING_TIMER(timer, "small sieve start positions");
                    /* When we're doing 2-level sieving, there is probably
                     * no real point in doing ssdpos initialization in
                     * several passes.
                     */
                    wss.ssd->small_sieve_prepare_many_start_positions(
                        first_region0_index,
                        std::min(SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE,
                                 ws.nb_buckets[1]),
                        ws.conf.logI, Q.sublat);
                    wss.ssd->small_sieve_activate_many_start_positions();
                },
                0);
        }

        pool.drain_queue(0);
        /* Now fill_in_buckets has completed for all levels. Time to check
         * that we had no overflow, and move on to process_bucket_region.
         */

        ws.check_buckets_max_full();
        auto exc = pool.get_exceptions<buckets_are_full>(0);
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
        // first_region0_index + ws.nb_buckets[ws.toplevel] * BRS[ws.toplevel] / BRS[1]);

        /* PROCESS THE REGIONS AT LEVEL 0 */
        process_many_bucket_regions(ws, wc_p, aux_p, Q, pool, first_region0_index, w);

        /* We need that, because the next downsort_tree call in the loop
         * above (for LEVEL>1) will reset the pointers while filling the 1l
         * buckets -- and we read the 1l buckets from PBR.
         */
        if (ws.toplevel > 1)
            pool.drain_queue(0);
    }
}

template <int LEVEL>
void downsort_tree(
    nfs_work & ws,
    std::shared_ptr<nfs_work_cofac> wc_p,
    std::shared_ptr<nfs_aux> aux_p,
    qlattice_basis const & Q,
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


#endif /* CADO_SIQS_FILL_IN_BUCKETS_INL */
