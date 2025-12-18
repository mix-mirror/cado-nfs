#ifndef CADO_SIQS_FILL_IN_BUCKETS_INL
#define CADO_SIQS_FILL_IN_BUCKETS_INL

void fill_in_buckets_prepare_plattices(
        MAYBE_UNUSED nfs_work & ws,
        MAYBE_UNUSED siqs_special_q_data const & Q,
        MAYBE_UNUSED thread_pool & pool,
        MAYBE_UNUSED int side,
        MAYBE_UNUSED cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS-1> & precomp_plattice)
{
}

template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static void fill_in_buckets_toplevel_sublat(
    bucket_array_t<LEVEL, TARGET_HINT> &,
    nfs_work &,
    siqs_special_q_data const &,
    fb_slice<FB_ENTRY_TYPE> const &,
    plattices_dense_vector_t *,
    where_am_I &)
{
    throw std::runtime_error("sublat is not supported in SIQS");
}

template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
void
fill_in_buckets_toplevel(bucket_array_t<LEVEL, TARGET_HINT> & orig_BA,
                         nfs_work & ws, fb_slice<FB_ENTRY_TYPE> const & slice,
                         siqs_special_q_data const & Q,
                         plattices_dense_vector_t * /* unused */,
                         where_am_I & w)
{
    int const logI = ws.conf.logI;
    size_t logJ = nbits(ws.J) - 1u; /* 2^m has m+1 bits */
    ASSERT_ALWAYS(ws.J == 1u << logJ); /* J must be a power of 2 */

    ASSERT_ALWAYS(!Q.sublat.m);

    /* local copy. Gain a register + use stack */
    bucket_array_t<LEVEL, TARGET_HINT> BA = std::move(orig_BA);

    slice_index_t const slice_index = slice.get_index();

    /* Write new set of pointers for the new slice */
    BA.add_slice_index(slice_index);
    WHERE_AM_I_UPDATE(w, i, slice_index);

    slice_offset_t hint = 0;

    int logB = LOG_BUCKET_REGIONS[LEVEL];
    typename bucket_array_t<LEVEL, TARGET_HINT>::update_t::br_index_t bmask =
        (1UL << logB) - 1;

    for (auto const & e: slice) {
        increment_counter_on_dtor<slice_offset_t> const _dummy(hint);
        if (!Q.is_coprime_to(e.p))
            continue;
        if (discard_power_for_bucket_sieving(e))
            continue;

        siqs_largesieve ple(Q, e, logJ);
        ASSERT_ALWAYS((ple.get_pp() >> logI) != 0); /* p must be > I */

        WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
        const fbprime_t p = slice.get_prime(hint);
        WHERE_AM_I_UPDATE(w, p, p);
#else
        const fbprime_t p = 0;
#endif

        typename bucket_array_t<LEVEL, TARGET_HINT>::update_t u(0, p, hint,
                                                                slice_index);

        std::vector<siqs_largesieve::T_elt> T1s;
        for (unsigned char i = 0; i < e.nr_roots; ++i) {
            ple.prepare_for_root(T1s, i, logI, 0u);
            auto const pp = ple.get_pp();
            auto it2 = ple.get_T2().rbegin();
            auto T2rend = ple.get_T2().rend(); /* T2 does not change during iteration */
            for (auto const & [t1, d1]: T1s) {

                /* case t1+t2 in [0, I[ */
                for (auto const & [t2, d2]: ple.get_T2()) {
                    if ((t1+t2) >> logI) {
                        break;
                    }
                    uint64_t x = ((d1 | d2) << logI) | (t1+t2);
                    u.set_x(x & bmask);
                    BA.push_update(x >> logB, u, w);
                }

                /* case t1+t2 in [p, p+I[ */
                if (it2 != T2rend && t1 + it2->first >= pp) {
                    while (it2 != T2rend && ((t1 + it2->first - pp) >> logI)) {
                        ++it2;
                    }
                    for (auto it = it2; it != T2rend; ++it) {
                        auto const [t2, d2] = *it;
                        if (t1+t2 < pp) {
                            break;
                        }
                        uint64_t x = ((d1 | d2) << logI) | (t1+t2-pp);
                        u.set_x(x & bmask);
                        BA.push_update(x >> logB, u, w);
                    }
                }
            }
        }
    }
    orig_BA = std::move(BA);
}

/* TARGET_HINT is shorthint_t or void */
template <int LEVEL, typename TARGET_HINT>
static void
fill_in_buckets_lowlevel(
        MAYBE_UNUSED bucket_array_t<LEVEL, TARGET_HINT> & orig_BA,
        MAYBE_UNUSED nfs_work & ws,
        MAYBE_UNUSED siqs_special_q_data const & Q,
        MAYBE_UNUSED plattices_vector_t & plattices_vector,
        MAYBE_UNUSED bool first_reg MAYBE_UNUSED,
        MAYBE_UNUSED where_am_I & w)
{
    throw std::runtime_error("LEVEL >= 2 sieving not implemented");
}

template <int LEVEL>
void downsort_tree(
        MAYBE_UNUSED nfs_work & ws,
        MAYBE_UNUSED std::shared_ptr<nfs_work_cofac> wc_p,
        MAYBE_UNUSED std::shared_ptr<nfs_aux> aux_p,
        MAYBE_UNUSED siqs_special_q_data const & Q,
        MAYBE_UNUSED thread_pool & pool,
        MAYBE_UNUSED uint32_t bucket_index, /* for the current level ! */
        MAYBE_UNUSED uint32_t first_region0_index,
        MAYBE_UNUSED std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> & precomp_plattices,
        MAYBE_UNUSED where_am_I & w)
{
    throw std::runtime_error("LEVEL >= 2 sieving not implemented");
}

#endif /* CADO_SIQS_FILL_IN_BUCKETS_INL */
