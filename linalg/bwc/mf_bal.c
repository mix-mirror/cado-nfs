#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <limits.h>
#include <sys/time.h>
#include "mf.h"
#include "mod_ul.h"
#include "balancing.h"
#include "rowset_heap.h"
#include "cheating_vec_init.h"
#include "mf_bal.h"
#include "portability.h" // strdup // IWYU pragma: keep
#include "timing.h"     // wct_seconds
#include "fix-endianness.h" // fread32_little
#include "macros.h"
#include "params.h"

typedef int (*sortfunc_t) (const void *, const void *);

int revcmp_2u32(const uint32_t * a, const uint32_t * b)
{
    if (a[0] > b[0]) return -1;
    if (b[0] > a[0]) return 1;
    if (a[1] > b[1]) return -1;
    if (b[1] > a[1]) return 1;
    return 0;
}

/* {{{ Datatype containing description for slices of the matrix */
struct slice {
    uint32_t * r;
    uint32_t nrows;
    uint32_t coeffs;
    uint32_t i0;
};

int slice_finding_helper(const uint32_t * key, const struct slice * elem)/*{{{*/
{
    if (*key < elem->i0) {
        return -1;
    } else if (*key >= elem->i0 + elem->nrows) {
        return 1;
    }
    return 0;
}

unsigned int which_slice(const struct slice * slices, unsigned int ns, uint32_t k)
{
    const struct slice * found = (const struct slice *) bsearch(
            (const void *) &k, (const void *) slices,
            ns, sizeof(struct slice),
            (sortfunc_t) &slice_finding_helper);
    if (found == NULL) {
        return UINT_MAX;
    } else {
        return found - slices;
    }
}/*}}}*/

void free_slices(struct slice * slices, unsigned int n)
{
    unsigned int i;
    for(i = 0 ; i < n ; i++) {
        free(slices[i].r);
    }
    free(slices);
}

struct slice * alloc_slices(unsigned int water, unsigned int n)
{
    struct slice * res;
    unsigned int i;

    res = (struct slice*) malloc(n * sizeof(struct slice));

    uint32_t common_size = water / n;

    // we now require equal-sized blocks.
    ASSERT_ALWAYS(water % n == 0);

    for(i = 0 ; i < n ; i++) {
        res[i].nrows = common_size; // min_size + (i < (water % n));
        res[i].r = (uint32_t *) malloc(res[i].nrows * sizeof(uint32_t));
        res[i].coeffs = 0;
    }

    return res;
}
/* }}} */
/* {{{ shuffle_rtable: This is the basic procedure which dispatches rows
 * (or columns) in buckets according to our preferred strategy
 */
struct slice * shuffle_rtable(
        const char * text,
        uint32_t (*rt)[2],
        uint32_t n,
        unsigned int ns)
{
    uint32_t i;
    struct bucket * heap;
    struct slice * slices;
    clock_t t;

    t = -clock();
    qsort(rt, n, 2 * sizeof(uint32_t), (sortfunc_t) &revcmp_2u32);
    t += clock();
    printf("sort time %.1f s\n", (double) t / CLOCKS_PER_SEC);

    
    t = -clock();

    slices = alloc_slices(n, ns);

    heap = (struct bucket *) malloc(ns * sizeof(struct bucket));
    for(i = 0 ; i < ns ; i++) {
        heap[i].s = 0;
        heap[i].i = i;
        heap[i].room = slices[i].nrows;
    }
    make_heap(heap, heap + ns);

    /* Then we're putting each row through the appropriate bucket, using
     * the heap to constantly update.
     * Eventually, we have constructed the set of rows going in the
     * bucket into slices[0]...slices[nslices-1]
     */
    for(i = 0 ; i < n ; i++) {
        pop_heap(heap, heap + ns);
        int j = heap[ns-1].i;
        int pos = slices[j].nrows-heap[ns-1].room;
        ASSERT(heap[ns-1].room);
        slices[j].r[pos] = rt[i][1];
        heap[ns-1].s += rt[i][0];
        heap[ns-1].room--;
        push_heap(heap, heap + ns);
    }

    t += clock();

    printf("heap fill time %.1f\n", (double) t / CLOCKS_PER_SEC);

    qsort(heap, ns, sizeof(struct bucket), (sortfunc_t) &heap_index_compare);

    for(i = 0 ; i < ns ; i++) {
        int j = heap[i].i;
        ASSERT(heap[i].i == (int) i);
        printf("%s slice %d, span=%ld, weight=%ld\n",
                text,
                i, slices[j].nrows - heap[i].room,
                heap[i].s);
        slices[j].coeffs = heap[i].s;
        if (heap[i].room != 0) {
            abort();
        }

#if 0
        // This is for giving the possibility of reading the matrix
        // linearly.
        if (sort_positionally) {
            qsort(slices[j].data, slices[j].nrows, sizeof(struct row),
                    (sortfunc_t) row_compare_index);
        }
#endif
    }
    free(heap);

    uint32_t i0 = 0;
    for(i = 0 ; i < ns ; i++) {
        slices[i].i0=i0;
        i0 += slices[i].nrows;
    }

    return slices;
}

/* }}} */

/* {{{ read the mfile header if we have it, deduce #coeffs */
void read_mfile_header(balancing_ptr bal, const char * mfile, int withcoeffs)
{
    struct stat sbuf_mat[1];
    int rc = stat(mfile, sbuf_mat);
    if (rc < 0) {
        fprintf(stderr, "Reading %s: %s (not fatal)\n", mfile, strerror(errno));
        printf("%s: %" PRIu32 " rows %" PRIu32 " cols\n",
                mfile, bal->h->nrows, bal->h->ncols);
        printf("%s: main input file not present locally, total weight unknown\n", mfile);
        bal->h->ncoeffs = 0;
    } else {
        bal->h->ncoeffs = sbuf_mat->st_size / sizeof(uint32_t) - bal->h->nrows;
        if (withcoeffs) {
            if (bal->h->ncoeffs & 1) {
                fprintf(stderr, "Matrix with coefficient must have an even number of 32-bit entries for all (col index, coeff). Here, %" PRIu64 " is odd.\n", bal->h->ncoeffs);
                abort();
            }
            bal->h->ncoeffs /= 2;
        }

        int extra = bal->h->ncols - bal->h->nrows;
        if (extra > 0) {
            printf( "%s: %" PRIu32 " rows %" PRIu32 " cols"
                    " (%d extra cols)"
                    " weight %" PRIu64 "\n",
                    mfile, bal->h->nrows, bal->h->ncols,
                    extra,
                    bal->h->ncoeffs);
        } else if (extra < 0) {
            printf( "%s: %" PRIu32 " rows %" PRIu32 " cols"
                    " (%d extra rows)"
                    " weight %" PRIu64 "\n",
                    mfile, bal->h->nrows, bal->h->ncols,
                    -extra,
                    bal->h->ncoeffs);
        } else {
            printf( "%s: %" PRIu32 " rows %" PRIu32 " cols"
                    " weight %" PRIu64 "\n",
                    mfile, bal->h->nrows, bal->h->ncols, bal->h->ncoeffs);
        }
    }
}
/* }}} */

void mf_bal_decl_usage(param_list_ptr pl)
{
   param_list_decl_usage(pl, "mfile", "matrix file (can also be given freeform)");
   param_list_decl_usage(pl, "rwfile", "row weight file (defaults to <mfile>.rw)");
   param_list_decl_usage(pl, "cwfile", "col weight file (defaults to <mfile>.cw)");
   param_list_decl_usage(pl, "out", "output file name (defaults to stdout)");
   param_list_decl_usage(pl, "quiet", "be quiet");
   param_list_decl_usage(pl, "rectangular", "accept rectangular matrices (for block Lanczos)");
   param_list_decl_usage(pl, "withcoeffs", "expect a matrix with explicit coefficients (not just 1s)");
   param_list_decl_usage(pl, "rowperm", "permute rows in priority (defaults to auto)");
   param_list_decl_usage(pl, "colperm", "permute columns in priority (defaults to auto)");
   param_list_decl_usage(pl, "skip_decorrelating_permutation", "solve for the matrix M instead of the matrix P*M with P a fixed stirring matrix");
}

void mf_bal_configure_switches(param_list_ptr pl, struct mf_bal_args * mba)
{
    param_list_configure_switch(pl, "--quiet", &mba->quiet);
    // param_list_configure_switch(pl, "--display-correlation", &display_correlation);
    param_list_configure_switch(pl, "--rectangular", &mba->rectangular);
    param_list_configure_switch(pl, "--withcoeffs", &mba->withcoeffs);
}


void mf_bal_parse_cmdline(struct mf_bal_args * mba, param_list_ptr pl, int * p_argc, char *** p_argv)
{
    unsigned int wild =  0;
    (*p_argv)++,(*p_argc)--;
    for(;(*p_argc);) {
        char * q;
        if (param_list_update_cmdline(pl, &(*p_argc), &(*p_argv))) continue;

        if ((*p_argv)[0][0] != '-' && wild == 0 && (q = strchr((*p_argv)[0],'x')) != NULL) {
            mba->nh = atoi((*p_argv)[0]);
            mba->nv = atoi(q+1);
            wild+=2;
            (*p_argv)++,(*p_argc)--;
            continue;
        }

        if ((*p_argv)[0][0] != '-' && wild == 0) { mba->nh = atoi((*p_argv)[0]); wild++,(*p_argv)++,(*p_argc)--; continue; }
        if ((*p_argv)[0][0] != '-' && wild == 1) { mba->nv = atoi((*p_argv)[0]); wild++,(*p_argv)++,(*p_argc)--; continue; }
        if ((*p_argv)[0][0] != '-' && wild == 2) {
            mba->mfile = (*p_argv)[0];
            wild++;
            (*p_argv)++,(*p_argc)--;
            continue;
        }
        fprintf(stderr, "unknown option %s\n", (*p_argv)[0]);
        exit(1);
    }
}

void mf_bal_interpret_parameters(struct mf_bal_args * mba, param_list_ptr pl)
{
    const char * tmp;

    if (!mba->nh || !mba->nv) {
        param_list_print_usage(pl, NULL, stderr);
        exit(EXIT_FAILURE);
    }

    if ((tmp = param_list_lookup_string(pl, "reorder")) != NULL) {
        if (strcmp(tmp, "auto") == 0) {
            mba->do_perm[1] = MF_BAL_PERM_AUTO;
            mba->do_perm[0] = MF_BAL_PERM_AUTO;
        } else if (strcmp(tmp, "rows") == 0) {
            mba->do_perm[1] = MF_BAL_PERM_NO;
            mba->do_perm[0] = MF_BAL_PERM_YES;
        } else if (strcmp(tmp, "columns") == 0) {
            mba->do_perm[1] = MF_BAL_PERM_YES;
            mba->do_perm[0] = MF_BAL_PERM_NO;
        } else if (strcmp(tmp, "rows,columns") == 0) {
            mba->do_perm[1] = MF_BAL_PERM_YES;
            mba->do_perm[0] = MF_BAL_PERM_YES;
        } else if (strcmp(tmp, "columns,rows") == 0) {
            mba->do_perm[1] = MF_BAL_PERM_YES;
            mba->do_perm[0] = MF_BAL_PERM_YES;
        } else if (strcmp(tmp, "both") == 0) {
            mba->do_perm[1] = MF_BAL_PERM_YES;
            mba->do_perm[0] = MF_BAL_PERM_YES;
        } else {
            fprintf(stderr, "Argument \"%s\" to the \"reorder\" parameter not understood\n"
                    "Supported values are:\n"
                    "\tauto (default)\n"
                    "\trows\n"
                    "\tcolumns\n"
                    "\tboth (equivalent forms: \"rows,columns\" or \"columns,rows\"\n",
                    tmp);
            exit(EXIT_FAILURE);
        }
    }

    if ((tmp = param_list_lookup_string(pl, "mfile")) != NULL) {
        mba->mfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "rwfile")) != NULL) {
        mba->rwfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "cwfile")) != NULL) {
        mba->cwfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "out")) != NULL) {
        mba->bfile = tmp;
    }
    param_list_parse_int(pl, "skip_decorrelating_permutation", &mba->skip_decorrelating_permutation);
}

void mf_bal_adjust_from_option_string(struct mf_bal_args * mba, const char * opts)
{
    /* Take a comma-separated list of balancing_options, and use that to
     * adjust the mba structure.
     */
    if (!opts) return;
    /* Create a new param_list from opts {{{ */
    char ** n_argv;
    char ** n_argv0;
    int n_argc;
    char * my_opts;
    ASSERT_ALWAYS(opts);
    my_opts = strdup(opts);
    n_argv0 = n_argv = malloc(strlen(my_opts) * sizeof(char*));
    n_argc = 0;
    n_argv[n_argc++]="mf_bal";
    for(char * q = my_opts, * qq; q != NULL; q = qq) {
        qq = strchr(q, ',');
        if (qq) { *qq++='\0'; }
        n_argv[n_argc++]=q;
    }

    param_list pl2;
    param_list_init(pl2);

    mf_bal_decl_usage(pl2);
    mf_bal_configure_switches(pl2, mba);
    mf_bal_parse_cmdline(mba, pl2, &n_argc, &n_argv);
    mf_bal_interpret_parameters(mba, pl2);

    param_list_clear(pl2);
    free(my_opts);
    free(n_argv0);

    /* }}} */
}

void mf_bal(struct mf_bal_args * mba)
{
    mf_bal_n(mba, 1);
}

void mf_bal_n(struct mf_bal_args * mbas, int n)
{
    const char * text[2] = { "row", "column", };
    int rc;
    int rectangular = mbas[0].rectangular;

    char *(*freeit)[2] = malloc(n * sizeof(char *[2]));
    balancing * bals = malloc(n * sizeof(balancing));

    /* {{{ get/build filenames */
    for(int midx = 0 ; midx < n ; midx++) {
        struct mf_bal_args * mba = &mbas[midx];

        /* we sometimes allocate strings, which need to be freed eventually */
        // char * freeit[2] = { NULL, NULL, };
        if (mba->mfile && !mba->rwfile) {
            mba->rwfile = freeit[midx][0] = build_mat_auxfile(mba->mfile, "rw", ".bin");
        }

        if (mba->mfile && !mba->cwfile) {
            mba->cwfile = freeit[midx][1] = build_mat_auxfile(mba->mfile, "cw", ".bin");
        }

        if (!mba->rwfile) { fprintf(stderr, "No rwfile given\n"); exit(1); }
        if (!mba->cwfile) { fprintf(stderr, "No cwfile given\n"); exit(1); }
        if (!mba->mfile) {
            fprintf(stderr, "Matrix file name (mfile) must be given, even though the file itself does not have to be present\n");
            exit(1);
        }
    }
    /*}}}*/

    /* {{{ init balancing structs */
    for(int midx = 0 ; midx < n ; midx++) {
        balancing_init(bals[midx]);
        bals[midx]->h->nh = mbas[midx].nh;
        bals[midx]->h->nv = mbas[midx].nv;
    }
    /* }}} */

    /* {{{ Get the row weights and column weights files (which *must exist*) */
    for(int midx = 0 ; midx < n ; midx++) {
        balancing_ptr bal = bals[midx];
        struct mf_bal_args * mba = &mbas[midx];

        struct stat sbuf[2][1];
        rc = stat(mba->rwfile, sbuf[0]);
        if (rc < 0) { perror(mba->rwfile); exit(1); }
        bal->h->nrows = sbuf[0]->st_size / sizeof(uint32_t);

        rc = stat(mba->cwfile, sbuf[1]);
        if (rc < 0) { perror(mba->cwfile); exit(1); }
        bal->h->ncols = sbuf[1]->st_size / sizeof(uint32_t);

        if (bal->h->ncols > bal->h->nrows) {
            fprintf(stderr, "Warning. More columns than rows. There could be bugs.\n");
        }

        read_mfile_header(bal, mba->mfile, mba->withcoeffs);
    }
    /*}}}*/

    /* {{{ Compute the de-correlating permutation.
     * This only makes sense on the last matrix (since it amounts to
     * changing it from M to M*S for some S)
     */
    {
        int midx = n - 1;
        balancing_ptr bal = bals[midx];
        struct mf_bal_args * mba = &mbas[midx];

        if (mba->skip_decorrelating_permutation) {
            /* internal, for debugging. This removes the de-correlating
             * permutation. Nothing to do with what is called
             * "shuffled-product" elsewhere, except that both are taken care
             * of within mf_bal. */
            bal->h->pshuf[0] = 1;
            bal->h->pshuf[1] = 0;
            bal->h->pshuf_inv[0] = 1;
            bal->h->pshuf_inv[1] = 0;
        } else {
            modulusul_t M;
            modul_initmod_ul(M, MIN(bal->h->nrows, bal->h->ncols));
            residueul_t a,b;
            residueul_t ai,bi;
            modul_init(a, M);
            modul_init(b, M);
            modul_init(ai, M);
            modul_init(bi, M);
            modul_set_ul(a, (unsigned long) sqrt(bal->h->ncols), M);
            modul_set_ul(b, 42, M);

            for( ; modul_inv(ai, a, M) == 0 ; modul_add_ul(a,a,1,M)) ;
            modul_mul(bi, ai, b, M);
            modul_neg(bi, bi, M);

            bal->h->pshuf[0] = modul_get_ul(a, M);
            bal->h->pshuf[1] = modul_get_ul(b, M);
            bal->h->pshuf_inv[0] = modul_get_ul(ai, M);
            bal->h->pshuf_inv[1] = modul_get_ul(bi, M);

            modul_clear(a, M);
            modul_clear(b, M);
            modul_clear(ai, M);
            modul_clear(bi, M);
            modul_clearmod(M);
        }
    }
    /* }}} */

    /* The grid size is given as (nh, nv).
     *
     * matrix 0 is split with grid (nh, nv)
     * matrix 1 is split with grid (nv, nh)
     * (and so on)
     *
     * so that if we have an odd number of matrices, we have the
     * (surprising) situation where
     * mbas[0].nh == mbas[n-1].nv
     *
     * (and in particular, the pair (mbas[0].nh, mbas[n-1].nv) is not
     * guaranteed to be (nh, nv))
     */
    unsigned int G = mbas[0].nh * mbas[0].nv;

    struct per_dim_stats {
        unsigned int matsize;
        unsigned int blocksize;
        unsigned int padding;
        uint32_t (*perm)[2];
        int todo;
    };

    struct per_matrix_stats {
        uint32_t * weights[2];
        uint32_t nzero[2];
        uint32_t totalweight;
        double avg[2];
        double sdev[2];
    };

    struct per_dim_stats * D = malloc((n + 1) * sizeof(struct per_dim_stats));
    struct per_matrix_stats * M = malloc(n * sizeof(struct per_matrix_stats));

    /* {{{ get (and check) the common dimensions */
    for(int midx = 0 ; midx < n ; midx++) {
        if ((midx & 1) == 0) {
            ASSERT_ALWAYS(mbas[midx].nh == mbas[0].nh);
            ASSERT_ALWAYS(mbas[midx].nv == mbas[0].nv);
        } else {
            ASSERT_ALWAYS(mbas[midx].nh == mbas[0].nv);
            ASSERT_ALWAYS(mbas[midx].nv == mbas[0].nh);
        }
        D[midx].matsize = bals[midx]->h->nrows;
        ASSERT_ALWAYS(!midx || bals[midx-1]->h->ncols == bals[midx]->h->nrows);
    }
    D[n].matsize = bals[n-1]->h->ncols;
    /* }}} */

    /* {{{ define block sizes and padding for both inner and outer dimensions */
    for(int didx = 0 ; didx < n + 1 ; didx++) {
        D[didx].blocksize = iceildiv(D[didx].matsize, G);

        /* We also want to enforce alignment of the block size with
         * respect to the SIMD things.
         * Given that mmt_vec_init provides 64-byte alignment of vector
         * areas, we may enforce the block size to be a multiple of
         * 8 in order to effectively guarantee 64-byte alignment for all
         * chunks. (Admittedly, this is a bit fragile; if we were to
         * possibly use smaller items, that would change stuff somewhat).
         */
        for ( ; D[didx].blocksize % (FORCED_ALIGNMENT_ON_MPFQ_VEC_TYPES / MINIMUM_ITEM_SIZE_OF_MPFQ_VEC_TYPES) ; D[didx].blocksize++);
    }

    if (!rectangular) {
        printf("Padding to a square matrix\n");
        D[0].blocksize = D[n].blocksize = MAX(D[0].blocksize, D[n].blocksize);
    }

    for(int didx = 0 ; didx < n + 1 ; didx++) {
        D[didx].padding = G * D[didx].blocksize - D[didx].matsize;
        printf("Padding dimension %d to %u+%u=%u , which is %u blocks of %u\n",
                didx, D[didx].matsize, D[didx].padding, G * D[didx].blocksize,
                G,
                D[didx].blocksize);
    }
    /* }}} */

    /* {{{ read weights files for all matrices, and do (and report) stats.
     * We also fill bal->h->{ncoeffs, nzrows, nzcols, flags} for each matrix
     */
    for(int midx = 0 ; midx < n ; midx++) {
        balancing_ptr bal = bals[midx];
        M[midx].totalweight = 0;
        for(int d = 0 ; d < 2 ; d++) {
            int didx = midx + d;
            size_t nitems = D[didx].matsize + D[didx].padding;

            M[midx].weights[d] = NULL;

            struct mf_bal_args * mba = &mbas[midx];
            const char * filename = d == 0 ? mba->rwfile : mba->cwfile;

            FILE * fw = fopen(filename, "rb");
            if (fw == NULL) { perror(filename); exit(1); }
            M[midx].weights[d] = malloc(nitems * sizeof(uint32_t));
            memset(M[midx].weights[d], 0, nitems * sizeof(uint32_t));
            /* Padding rows and cols have zero weight of course */
            double t_w;
            t_w = -wct_seconds();
            size_t nr = fread32_little(M[midx].weights[d], D[didx].matsize, fw);
            t_w += wct_seconds();
            fclose(fw);
            
            if (nr < D[didx].matsize) {
                fprintf(stderr, "%s: short %s count\n", filename, text[d]);
                exit(1);
            }
            printf("read %s in %.1f s (%.1f MB / s)\n", filename, t_w,
                    1.0e-6 * nr * sizeof(uint32_t) / t_w);

            /* {{{ compute and report average weight and sdev */
            {
                double s1 = 0;
                double s2 = 0;
                uint64_t totalweight = 0;
                for(unsigned int r = 0 ; r < D[didx].matsize ; r++) {
                    double x = M[midx].weights[d][r];
                    totalweight += M[midx].weights[d][r];
                    s1 += x;
                    s2 += x * x;
                }
                /* {{{ make sure we agree on the total weight */
                if (M[midx].totalweight) {
                    /* hopefully the sums of row weights and column
                     * weights coincide!
                     */
                    ASSERT_ALWAYS(M[midx].totalweight == totalweight);
                } else {
                    M[midx].totalweight = totalweight;
                }
                if (bal->h->ncoeffs) {
                    if (totalweight != bal->h->ncoeffs) {
                        fprintf(stderr,
                                "Inconsistency in number of coefficients\n"
                                "From %s: %" PRIu64
                                ", from file sizes; %" PRIu64 "\n",
                                filename, totalweight, bal->h->ncoeffs);
                        fprintf(stderr,
                                "Maybe use the --withcoeffs option for DL matrices ?\n");
                        exit(EXIT_FAILURE);
                    }
                } else {
                    bal->h->ncoeffs = totalweight;
                    printf("%" PRIu64 " coefficients counted\n", totalweight);
                }
                /*}}}*/
                M[midx].avg[d] = s1 / D[didx].matsize;
                M[midx].sdev[d] = sqrt(s2 / D[didx].matsize - M[midx].avg[d]*M[midx].avg[d]);
                printf("%" PRIu32 " %ss ;"
                        " avg %.1f sdev %.1f"
                        " [scan time %.1f s]\n",
                        D[didx].matsize, text[d],
                        M[midx].avg[d], M[midx].sdev[d],
                        t_w);
            }
            /* }}} */

            /* {{{ count zero rows or columns */
            M[midx].nzero[d] = 0;
            for(size_t r = 0 ; r < D[didx].matsize ; r++) {
                M[midx].nzero[d] += (M[midx].weights[d][r] == 0);
            }
            /* }}} */

        }
        bal->h->nzrows = M[midx].nzero[0];
        bal->h->nzcols = M[midx].nzero[1];
        bal->h->flags = 0;
    }
    /* }}} */

    /* {{{ Determine what we have to do with each dimension. */
    /* Each dimension exists for two matrices, except for the
     * outer ones which of course exist only for one.
     * At times, for a given index didx inside the list of n+1 matrix
     * dimensions ([0] to [n]), we need to consider the following
     * matrices:
     *  - the "previous" ("left") matrix is for d=1, and it's matrix
     *  number didx-1 if didx>0. (we're interested in its columns,
     *  hence d=1)
     *  - the "next" ("right") matrix is for d=0, and it's matrix
     *  number didx, if didx<n (that is, unless we're speaking of the last
     *  dimension, which is numbered [n] and corresponds to the
     *  number of columns of matrix [n-1])
     */

    /* For outer dimensions, we want to compute a balancing permutation
     * based on the single set of weights we have access to.
     *
     * For inner dimensions (inner dimensions happen only when we have a
     * chain of several matrices), there's a decision to make about what
     * is the "weight" we should rely on: we have two !
     *
     * Short of something smarter, a simple proposal is to do exactly as
     * we do with the outer dimensions, in fact: we look at the two
     * standard deviations, and permute according to whichever is larger.
     */
    for(int didx = 0 ; didx < n + 1 ; didx++) {
        int dp[2];
        int mi[2];
        for(int d = 0 ; d < 2 ; d++)  {
            mi[d] = (didx + n - d) % n;
            dp[d] = mbas[mi[d]].do_perm[d];
        }
        if (dp[0] == MF_BAL_PERM_YES && dp[1] == MF_BAL_PERM_NO) {
            /* ok */
        } else if (dp[0] == MF_BAL_PERM_NO && dp[1] == MF_BAL_PERM_YES) {
            /* ok */
        } else if (dp[0] == MF_BAL_PERM_NO && dp[1] == MF_BAL_PERM_NO) {
            /* ok. No work to do here ! */
        } else if (dp[0] == MF_BAL_PERM_AUTO && dp[1] == MF_BAL_PERM_AUTO) {
            int choose = M[mi[1]].sdev[1] > M[mi[0]].sdev[0];
            printf("Choosing a %s (matrix %d) permutation based"
                    " on largest deviation"
                    " (%.2f > %.2f)\n",
                    text[choose], mi[choose],
                    M[mi[choose]].sdev[choose],
                    M[mi[!choose]].sdev[!choose]);
            mbas[mi[choose]].do_perm[choose] = MF_BAL_PERM_YES;
            mbas[mi[!choose]].do_perm[!choose] = MF_BAL_PERM_NO;
        } else {
            /* in particular, AUTO only works with AUTO, and YES
             * obviously does not work with YES.
             */
            fprintf(stderr,
                    "inconsistent requests for balancing order"
                    " at dimension %d ; doperm[%d][1]=%d, doperm[%d][0]=%d\n",
                    didx,
                    didx-1, mbas[mi[1]].do_perm[1],
                    didx, mbas[mi[0]].do_perm[0]);
            exit(EXIT_FAILURE);
        }
        /* At this point we only have YESes and Nos */
    }
    /* }}} */


    /* {{{ Prepare the tables that we're going to sort */
    for(int didx = 0 ; didx < n + 1 ; didx++) {
        size_t nitems = D[didx].matsize + D[didx].padding;

        D[didx].perm = NULL;

        /* Are we going to permute this dimension based on the rows of
         * the next matrix, or the columns of the previous one ?
         */
        /* We let didx run from 0 to n included, but in reality if the
         * product of the whole matrix chain is square, dimensions [0]
         * and [n] are the same, and in that case we're comparing the
         * same things.
         *
         * I.e., didx=0 compares rows of M (d=0) with columns of M (d=1)
         *   and didx=1 also compares rows of M (d=0) with columns of M (d=1)
         */
        uint32_t * active_weights = NULL;
        for(int d = 0 ; d < 2 ; d++)  {
            int midx = (didx + n - d) % n;
            if (mbas[midx].do_perm[d] == MF_BAL_PERM_YES) {
                ASSERT_ALWAYS(active_weights == NULL);
                active_weights = M[midx].weights[d];
            }
        }

        if (!active_weights) {
            /* we're not rearranging this table */
            continue;
        }

        D[didx].perm = malloc(nitems * sizeof(uint32_t[2]));

        /* prepare for qsort */
        for(size_t r = 0 ; r < nitems ; r++) {
            /* Compute rx so that the column r in the matrix we work with
             * is actually column rx in the original matrix.  */

            /* We do so a priori only for the last dimension only */
            int shuffle = didx == n;

            /* so the base case is simply rx == r */
            size_t rx = r;

#if 1   /* very weird behaviour, to be investigated */
            /* FIXME: previous code did that both for d==1 *AND* d==0,
             * and I very much think that it's bogus to do that for d==0.
             * It could be that it's a mistake that has gone unnoticed.
             */
            if (n == 1)
                shuffle = 1;
#endif

            if (shuffle && r < D[didx].matsize) {
                /* The balancing_pre_* function satisfy
                 * f([0,bal->h->ncols[) \subset [0,bal->h->ncols[.
                 * and correspond to identity for x>=bal->h->ncols
                 */
                rx = balancing_pre_unshuffle(bals[n-1], r);
                ASSERT(balancing_pre_shuffle(bals[n-1], rx) == r);
                ASSERT(rx < D[didx].matsize);
            }
            D[didx].perm[r][0] = active_weights[rx];
            D[didx].perm[r][1] = r;
        }
    }
    /* }}} */

#if 0 /* {{{ examine row/col correlation (disabled) */
    /* We used to examine the correlation between row and column
     * weight. This is important, as it accounts for some timing jitter
     * on the local matrix products.
     *
     * Unfortunately, in full generality, I'm having difficulties to make
     * sense out of this notion for rectangular matrices, so let's
     * comment it out for the moment.
     */

    if (display_correlation) {
        size_t nitems = bal->h->nrows + rpadding;
        FILE * frw = fopen(rwfile, "rb");
        if (frw == NULL) { perror(rwfile); exit(1); }
        uint32_t * rowweights = malloc(nitems * sizeof(uint32_t));
        memset(rowweights, 0, nitems * sizeof(uint32_t));
        double t_rw;
        t_rw = -wct_seconds();
        size_t nr = fread32_little(rowweights, bal->h->nrows, frw);
        t_rw += wct_seconds();
        fclose(frw);
        if (nr < bal->h->nrows) {
            fprintf(stderr, "%s: short row count\n", rwfile);
            exit(1);
        }
        double rs1 = 0;
        double rs2 = 0;
        double rc_plain = 0;
        double rc_decorr = 0;
        for(size_t r = 0 ; r < nitems ; r++) {
            double x = rowM[r].weights;
            rs1 += x;
            rs2 += x * x;
            rc_plain += x * colM[r].weights;
            rc_decorr += x * bal->colperm[2*r];
        }
        double ravg = rs1 / nitems;
        double rsdev = sqrt(rs2 / nitems - ravg*ravg);
        double cavg = s1 / nitems;
        double csdev = sqrt(s2 / nitems - cavg*cavg);
        double pcov = rc_plain / nitems - ravg * cavg;
        double pcorr = pcov / csdev / rsdev;
        double dcov = rc_decorr / nitems - ravg * cavg;
        double dcorr = dcov / csdev / rsdev;

        printf("%" PRIu32 " rows ; avg %.1f sdev %.1f [scan time %.1f s]\n",
                bal->h->nrows, ravg, rsdev, t_rw);
        printf("row-column correlation coefficient is %.4f\n",
                pcorr);
        printf("row-column correlation coefficient after decorrelation is %.4f\n",
                dcorr);
    }
#endif/*}}}*/

    /* {{{ All weight tables can be freed now */
    for(int midx = 0 ; midx < n ; midx++) {
        for(int d = 0 ; d < 2 ; d++) {
            free (M[midx].weights[d]);
        }
    }
    /* }}} */

    for(int didx = 0 ; didx < n + 1 ; didx++) {
        size_t nitems = D[didx].matsize + D[didx].padding;
        for(int d = 0 ; d < 2 ; d++)  {
            int midx = (didx + n - d) % n;

            if (mbas[midx].do_perm[d] == MF_BAL_PERM_NO)
                /* in effect, only one value of d will run through the
                 * loop */
                continue;

            balancing_ptr bal = bals[midx];
            bal->h->flags |= d == 0 ? FLAG_ROWPERM : FLAG_COLPERM;
            uint32_t **pperm = d == 0 ? &bal->rowperm : &bal->colperm;

            // unsigned int k = G; // better choice, perhaps ?
            unsigned int k = d == 0 ? mbas[midx].nh : mbas[midx].nv;

            struct slice * h = shuffle_rtable(text[d], D[didx].perm, nitems, k);
            *pperm = malloc(nitems * sizeof(uint32_t));

            for(unsigned int ii = 0 ; ii < k ; ii++) {
                const struct slice * r = &(h[ii]);
                memcpy(*pperm + r->i0, r->r, r->nrows * sizeof(uint32_t));
            }
            free_slices(h, k);
        }
    }

    /*
     * "replicating" a permutation means connecting the rows
     * of the first matrix with the columns of the last one.
     *
     * In the multi matrix case, we're not quite sure that we want to use
     * this flag, but for consistency with the behaviour of the single
     * matrix code, we keep the old behaviour that if the matrix is
     * square, we forcibly do this replication
     */
    if (n == 1 && !rectangular) {
        bals[0]->h->flags |= FLAG_REPLICATE;
    }

    for(int midx = 0 ; midx < n ; midx++) {
        balancing_ptr bal = bals[midx];
        struct mf_bal_args * mba = &mbas[midx];

        balancing_finalize(bal);

        balancing_write(bal, mba->mfile, mba->bfile);
        balancing_clear(bal);

        if (freeit[midx][0]) free(freeit[midx][0]);
        if (freeit[midx][1]) free(freeit[midx][1]);
    }
    free(freeit);
    free(bals);
    free(D);
    free(M);
}

