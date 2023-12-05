/* merge --- new merge program, for double-matrix algorithm

TODO:
* fix the discrepancy between the weight of R printed by merge (R has W=...)
  and that printed by replay-dblemat (Weight of the sparse submatrix: ...)
* fix the total weight W in merge, which should correspond to what replay says
  (ok for the c60 from README, and for
  parameters/polynomials/{c70,c80,c90,c100}.poly,
  but we get a discrepancy for the c110 from
  parameters/polynomials/c110.poly, where on baguette.loria.fr merge gives
  W=48426850 but replay gives 48345848).

Copyright 2019-2023 Charles Bouillaguet and Paul Zimmermann.

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* This code implements the algorithm described in reference [1].
   Reference [2] might be useful too.

   [1] Parallel Structured Gaussian Elimination for the Number Field Sieve,
       Charles Bouillaguet and Paul Zimmermann, Mathematical Cryptology,
       volume 0, number 1, pages 22-39, 2020.
   [2] Design and Implementation of a Parallel Markowitz Threshold Algorithm,
       Timothy A. Davis, Iain S. Duff, and Stojce Nakov
       SIAM Journal on Matrix Analysis and Applications
       Volume 41, Issue 2, 2020, https://doi.org/10.1137/19M1245815.

   The output format ("history") is described in README.fileformat
*/

#include "cado.h" // IWYU pragma: keep
/* the following should come after cado.h, which sets -Werror=all */
#ifdef  __GNUC__
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>       // for PRIu64
#include <stdint.h>         // for uint64_t, int32_t, int64_t, uint32_t
#include <string.h>         // for memset, memcpy, strlen
#ifdef HAVE_MINGW
#include <fcntl.h>         /* for _O_BINARY */
#endif
#include "filter_io.h"  // earlyparsed_relation_ptr
#ifdef FOR_DL
#include "gcd.h"
#endif
#include "gzip.h"       // fopen_maybe_compressed
#include "macros.h"
#include "memory.h"             // malloc_aligned
#include "memusage.h"   // PeakMemusage
#include "misc.h"       // UMAX
#include "mst.h"
#include "omp_proxy.h"
#include "params.h"     // param_list_parse_*
#include "purgedfile.h"     // for purgedfile_read_firstline
#include "sparse.h"
#include "timing.h"  // seconds
#include "typedefs.h"  // weight_t
#include "verbose.h"    // verbose_interpret_parameters
#include "merge_heap.h"
#include "merge_bookkeeping.h"
#include "merge_compute_weights.h"
#include "read_purgedfile_in_parallel.h"

/* MERGE_STRATEGY=1: classical strategy, where we optimize the weight increase
   in the matrix product M = L*R
   MERGE_STRATEGY=2: alternate strategy, where we optimize the weight increase
   in both L and R */
#define MERGE_STRATEGY 1

#if !(MERGE_STRATEGY == 1 || MERGE_STRATEGY == 2)
#error "MERGE_STRATEGY should be 1 or 2"
#endif

int twomerge_mode = 1;   // stays true until we do the first k-merge with k > 2

#ifdef DEBUG
static void
fprint_row (FILE *fp, filter_matrix_t *mat, index_t i)
{
  ASSERT_ALWAYS(mat->rows[i] != NULL);
  fprintf (fp, "%u", matLengthRow(mat, i));
  for (index_t k = 1; k <= matLengthRow(mat, i); k++)
    fprintf (fp, " %u", rowCell(mat->rows[i],k));
  fprintf (fp, "\n");
}
#endif


/*************************** output buffer ***********************************
 * There is one buffer per thread, so they can be accessed without synchronization
 * They are flushed in-order, so it is important that operations logged during a
 * single pass commute. In fact they do because the eliminated columns do not share a row
 */ 

typedef struct {
  char* buf;
  size_t size;  /* used size */
  size_t alloc; /* allocated size */
} buffer_struct_t;

static buffer_struct_t*
buffer_init (int nthreads)
{
  buffer_struct_t *Buf;
  Buf = malloc (nthreads * sizeof (buffer_struct_t));
  for (int i = 0; i < nthreads; i++) {
      Buf[i].buf = NULL;
      Buf[i].size = 0;
      Buf[i].alloc = 0;
    }
  return Buf;
}

static void
buffer_add (buffer_struct_t *buf, char *s)
{
  size_t n = strlen (s) + 1; /* count final '\0' */
  if (buf->size + n > buf->alloc) {
      buf->alloc = buf->size + n;
      buf->alloc += buf->alloc / MARGIN;
      buf->buf = realloc (buf->buf, buf->alloc * sizeof (char));
    }
  memcpy (buf->buf + buf->size, s, n * sizeof (char));
  buf->size += n - 1; /* don't count final '\0' */
}

static void
buffer_flush (buffer_struct_t *Buf, int nthreads, FILE *out)
{
  double cpu = seconds (), wct = wct_seconds ();
  for (int i = 0; i < nthreads; i++)
    {
      /* it is important to check whether size=0, otherwise the previous
         buffer will be printed twice */
      if (Buf[i].size != 0)
        fprintf (out, "%s", Buf[i].buf);
      Buf[i].size = 0;
    }
  cpu = seconds () - cpu;
  wct = wct_seconds () - wct;
  print_timings ("   buffer_flush took", cpu, wct);
  cpu_t[FLUSH] += cpu;
  wct_t[FLUSH] += wct;
}

static void
buffer_clear (buffer_struct_t *Buf, int nthreads)
{
  for (int i = 0; i < nthreads; i++)
    free (Buf[i].buf);
  free (Buf);
}

/*****************************************************************************/

static void
declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "mat", "input purged file");
  param_list_decl_usage(pl, "out", "output history file");
  param_list_decl_usage(pl, "skip", "number of heavy columns to bury (default "
				    CADO_STRINGIZE(DEFAULT_MERGE_SKIP) ")");
  param_list_decl_usage(pl, "target_density", "stop when the average row density exceeds this value"
			    " (default " CADO_STRINGIZE(DEFAULT_MERGE_TARGET_DENSITY) ")");
  param_list_decl_usage(pl, "force-posix-threads", "force the use of posix threads, do not rely on platform memory semantics");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
  param_list_decl_usage(pl, "t", "number of threads");
  param_list_decl_usage(pl, "v", "verbose mode");
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

/* check that mat->tot_weight and mat->wt say the same thing.
 *
 * Note that this check makes sense only if col_weight_t is a wide enough
 * type.  Either we have col_weight_t very inaccurate, and then the
 * invariant can't be right, or we choose to use a wider type for
 * col_weight_t, which should go with more accurate tracking in mat->wt,
 * and eventually the invariant should be always right
 */
void check_invariant(filter_matrix_t *mat)
{
    if (sizeof(col_weight_t) == 1)
      return;

    uint64_t tot_weight2 = 0;
    for (index_t i = 0; i < mat->ncols; i++) {
        ASSERT_ALWAYS(mat->wt[i] <= mat->rem_nrows);
        tot_weight2 += mat->wt[i];
    }
    printf("invariant %s: tw = %" PRIu64 " tw2=%" PRIu64 "\n",
            ok_NOK(mat->tot_weight == tot_weight2),
            mat->tot_weight, tot_weight2);
    
    ASSERT_ALWAYS(mat->tot_weight == tot_weight2);
}

#ifndef FOR_DL
/* sort row[0], row[1], ..., row[n-1] in non-decreasing order */
static void
sort_relation (index_t *row, unsigned int n)
{
  unsigned int i, j;

  for (i = 1; i < n; i++)
    {
      index_t t = row[i];
      if (t < row[i-1])
	{
	  row[i] = row[i-1];
	  for (j = i - 1; j > 0 && t < row[j-1]; j--)
	    row[j] = row[j-1];
	  row[j] = t;
	}
    }
}
#endif

/* callback function called by filter_rels */
static void *
insert_rel_into_table (void *context_data, earlyparsed_relation_ptr rel)
{
  filter_matrix_t *mat = (filter_matrix_t *) context_data;
  unsigned int j = 0;
  typerow_t buf[UMAX(weight_t)]; /* rel->nb is of type weight_t */

  for (unsigned int i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
    /* we no longer touch mat->wt, mat->rem_ncols, and mat->tot_weight
     * from here ; see compute_weights
     */
    if (h < mat->skip)
	continue; /* we skip (bury) the first 'skip' indices */
#ifdef FOR_DL
    exponent_t e = rel->primes[i].e;
    /* For factorization, they should not be any multiplicity here.
       For DL we do not want to count multiplicity in mat->wt */
    buf[++j] = (ideal_merge_t) {.id = h, .e = e};
#else
    ASSERT(rel->primes[i].e == 1);
    buf[++j] = h;
#endif
  }

#ifdef FOR_DL
  buf[0].id = j;
#else
  buf[0] = j;
#endif

  /* sort indices to ease row merges */
#ifndef FOR_DL
  sort_relation (&(buf[1]), j);
#else
  qsort (&(buf[1]), j, sizeof(typerow_t), cmp_typerow_t);
#endif

  mat->rows[rel->num] = heap_alloc_row(mat->heap, rel->num, j);
  compressRow (mat->rows[rel->num], buf, j);  /* sparse.c, simple copy loop... */

  return NULL;
}


static void
filter_matrix_read (filter_matrix_t *mat, const char *purgedname)
{
  uint64_t nread;
  char *fic[2] = {(char *) purgedname, NULL};

  /* first check if purgedname is seekable. if yes, we can do multithread
   * I/O */
  int can_go_parallel;
  {
      FILE * f = fopen_maybe_compressed(purgedname, "r");
      ASSERT_ALWAYS(f != NULL);
      can_go_parallel = fseek(f, 0, SEEK_END) == 0;
      fclose_maybe_compressed (f, purgedname);
  }

  if (!can_go_parallel) {
      fprintf(stderr, "# cannot seek in %s, using single-thread I/O\n", purgedname);
      /* read all rels */
      nread = filter_rels (fic, (filter_rels_callback_t) &insert_rel_into_table,
              mat, EARLYPARSE_NEED_INDEX_SORTED, NULL, NULL);
  } else {
      nread = read_purgedfile_in_parallel(mat, purgedname);
  }

  ASSERT_ALWAYS(nread == mat->nrows);
  mat->rem_nrows = nread;
}

/* check the matrix rows are sorted by increasing index */
/* this should not be needed at all, now that we ask filter_rels to
 * give us sorted rows with EARLYPARSE_NEED_INDEX_SORTED. (non-sorted
 * rows are sorted on the fly if needed)
 */
static void
check_matrix (filter_matrix_t *mat)
{
  #pragma omp parallel for
  for (index_t i = 0; i < mat->nrows; i++)
    {
      index_t l = matLengthRow (mat, i);
      for (index_t j = 1; j < l; j++)
        /* for DL we can have duplicate entries in the purged file, but after
           filter_matrix_read() they should be accumulated into one single
           entry (k,e), thus successive values of k cannot be equal here */
        if (matCell (mat, i, j) >= matCell (mat, i, j+1))
          {
            fprintf (stderr, "Error, the rows of the purged file should be sorted by increasing index\n");
            exit (EXIT_FAILURE);
          }
    }
}

/* Put the initial weights (after the initial 2-merges) in R->wt[].
   Warning: these weights are saturated at 255, thus cannot be used
   to compute the total weight of R. */
static void
compute_weights_R (filter_matrix_t *R, filter_matrix_t *mat)
{
#pragma omp parallel for
  for (index_t j = 0; j < mat->ncols; j++)
    R->wt[j] = mat->wt[j];
}

/* Stack non-empty columns at the beginning. Update mat->p (for DL) and jmin.
   We also have to "recompress" R->wt[]. */
static void recompress(filter_matrix_t *R, filter_matrix_t *mat, index_t *jmin)
{
	double cpu = seconds (), wct = wct_seconds ();
	uint64_t nrows = mat->nrows;
	uint64_t ncols = mat->ncols;

	/* sends the old column number to the new one */
	index_t *p = malloc(ncols * sizeof(*p));

        /* new column weights */
        col_weight_t *nwt = malloc(mat->rem_ncols * sizeof(*nwt));
        col_weight_t *nwtr = malloc(mat->rem_ncols * sizeof(*nwt));

        /* compute the number of non-empty columns */
        {
            /* need an array that is visible to all threads in order to do the
             * prefix sum. Wish I knew another way.
             */
            int T = omp_get_max_threads();
            index_t tm[T]; /* #non-empty columns seen by thread t */
#pragma omp parallel
            {
                int T = omp_get_num_threads();
                int t = omp_get_thread_num();
                index_t m = 0;
#pragma omp for schedule(static) nowait /* static is mandatory here */
                for (index_t j = 0; j < ncols; j++)
                    if (0 < mat->wt[j])
                        m++;
                tm[t] = m;

#pragma omp barrier

                /* prefix-sum over the T threads (sequentially) */
#pragma omp single
                {
                    index_t s = 0;
                    for (int t = 0; t < T; t++) {
                        index_t m = tm[t];
                        tm[t] = s;
                        s += m;
                    }
                    /* we should have s = mat->rem_ncols now, thus no need
                       to copy s into mat->rem_ncols, but it appears in
                       some cases it does not hold (cf
https://cado-nfs-ci.loria.fr/ci/job/future-parallel-merge/job/compile-debian-testing-amd64-large-pr/147/) */
                    mat->rem_ncols = s;
                }

                /* compute the new column indices */
                m = tm[t];
#pragma omp for schedule(static) /* static is mandatory here */
                for (index_t j = 0; j < ncols; j++) {
                    ASSERT(m <= j);
                    p[j] = m;
                    if (0 < mat->wt[j])
                        m++;
                }

                /* rewrite the row indices */
#pragma omp for schedule(guided) /* guided is slightly better than static */
                for (uint64_t i = 0; i < nrows; i++) {
                    if (mat->rows[i] == NULL) 	/* row was discarded */
                        continue;
                    for (index_t l = 1; l <= matLengthRow(mat, i); l++)
                        matCell(mat, i, l) = p[matCell(mat, i, l)];
                }

                /* update mat->wt and R->wt */
#pragma omp for schedule(static) /* static is slightly better than guided */
                for (index_t j = 0; j < ncols; j++)
                    if (0 < mat->wt[j])
                    {
                      nwt[p[j]] = mat->wt[j];
                      /* Note: during the 2-merges, R->wt[] is not
                         initialized, thus the following is useless
                         (but does not hurt).
                         After the 2-merges, we call
                         compute_weights_R() which initializes R->wt[]. */
                      nwtr[p[j]] = R->wt[j];
                    }

            } /* end parallel section */
        }

        /* We keep the inverse of p, to print the
	   original columns in the history file.
	   Warning: for a column j of weight 0, we have p[j] = p[j'] where
	   j' is the smallest column > j of positive weight, thus we only consider
	   j such that p[j] < p[j+1], or j = ncols-1. */
        if (mat->p == NULL) {
        	mat->p = malloc(mat->rem_ncols * sizeof (index_t));
                /* We must pay attention to the case of empty columns at the
                 * end */
		for (uint64_t i = 0, j = 0; j < mat->ncols && i < mat->rem_ncols; j++) {
                    if (p[j] == i && (j + 1 == mat->ncols || p[j] < p[j+1]))
                        mat->p[i++] = j; /* necessarily i <= j */
                }
        } else {
	/* update mat->p. It sends actual indices in mat to original indices in the purge file */
        // before : mat->p[i] == original
        //  after : mat->p[p[i]] == original
        /* Warning: in multi-thread mode, one should take care not to write
           some mat->p[j] before it is used by another thread.
           Consider for example ncols = 4 with 2 threads, and active
           columns 1 and 2. Then we have p[1] = 0 and p[2] = 1.
           Thus thread 0 executes mat->p[0] = mat->p[1], and thread 1 executes
           mat->p[1] = mat->p[2]. If thread 1 is ahead of thread 0, the final
           value of mat->p[0] will be wrong (it will be the initial value of
           mat->p[2], instead of the initial value of mat->p[1]).
           To solve that problem, we store the new values in another array. */
        	index_t *new_p = malloc (mat->rem_ncols * sizeof (index_t));
		/* static slightly better than guided for the following loop */
                #pragma omp for schedule(static)
		for (index_t j = 0; j < ncols; j++)
			if (0 < mat->wt[j])
				new_p[p[j]] = mat->p[j];
		free(mat->p);
		mat->p = new_p;
        }
        
	free(mat->wt);
	mat->wt = nwt;
        free(R->wt);
        R->wt = nwtr;

	/* update jmin */
	if (jmin[0] == 1)
                for (int w = 1; w <= MERGE_LEVEL_MAX; w++)
                /* Warning: we might have jmin[w] = ncols. */
                        jmin[w] = (jmin[w] < ncols) ? p[jmin[w]] : mat->rem_ncols;

	free(p);

	/* this was the goal all along! */
	mat->ncols = mat->rem_ncols;
	cpu = seconds () - cpu;
	wct = wct_seconds () - wct;
	print_timings ("   recompress took", cpu, wct);
	cpu_t[RECOMPRESS] += cpu;
	wct_t[RECOMPRESS] += wct;
}


/* For 1 <= w <= MERGE_LEVEL_MAX, put in jmin[w] the smallest index j such that
   mat->wt[j] = w. This routine is called only once, at the first call of
   compute_weights. */

static void
compute_jmin (filter_matrix_t *mat, index_t *jmin)
{
    {
        /* unfortunately, reduction on array sections requires OpenMP >= 4.5,
           which is not yet THAT widespread. We work around the problem */
        /* TODO: I wonder which openmp level we require anyway. Maybe
         * it's already 4.5+ */
        index_t tjmin[omp_get_max_threads()][MERGE_LEVEL_MAX + 1];

#pragma omp parallel /* reduction(min: jmin[1:MERGE_LEVEL_MAX]) */
        {
            int T = omp_get_num_threads();
            int tid = omp_get_thread_num();

            index_t *local = tjmin[tid];

            /* first initialize to ncols */
            for (int w = 1; w <= MERGE_LEVEL_MAX; w++)
                local[w] = mat->ncols;

	    /* compute_jmin takes so little time that it makes no sense
	       optimizing the schedule below */
            #pragma omp for schedule(static)
            for (index_t j = 0; j < mat->ncols; j++) {
                col_weight_t w = mat->wt[j];
                if (0 < w && w <= MERGE_LEVEL_MAX && j < local[w])
                    local[w] = j;
            }

	    /* compute_jmin takes so little time that it makes no sense
	       optimizing the schedule below */
            #pragma omp for schedule(static)
            for (int w = 1; w <= MERGE_LEVEL_MAX; w++) {
                jmin[w] = mat->ncols;
                for (int t = 0; t < T; t++)
                    if (jmin[w] > tjmin[t][w])
                        jmin[w] = tjmin[t][w];
            }
        }
    }

  jmin[0] = 1; /* to tell that jmin was initialized */

  /* make jmin[w] = min(jmin[w'], 1 <= w' <= w) */
  for (int w = 2; w <= MERGE_LEVEL_MAX; w++)
    if (jmin[w - 1] < jmin[w])
      jmin[w] = jmin[w - 1];
}

/* 
 * This does a pass on the matrix data (all rows), and collects the
 * following info
 *
 * mat->wt[]
 * mat->rem_ncols
 * mat->tot_weight
 *
 * In the general case, column weights in mat->wt need only be computed
 * up to cwmax + 1 since we only need to know whether the weights are <=
 * cwmax or not).
 *
 * However, this is not true for the shrink case, where a full count is
 * needed in order to accurately compute the density.
 *
 */
static void
compute_weights (filter_matrix_t *mat, index_t *jmin)
{
  double cpu = seconds (), wct = wct_seconds ();
  // col_weight_t cwmax = mat->cwmax;

  /* This function used to work with jmin already initialized, and maybe
   * still does. The thing is that it hasn't been used this way for a
   * while, and the call path with jmin[0] != 0 is not tested at all. If
   * needed, remove this assert, but be cautious !
   */
  ASSERT_ALWAYS(jmin[0] == 0);

  index_t j0;
  if (jmin[0] == 0) /* jmin was not initialized */
    {
      j0 = 0;
      // cwmax = MERGE_LEVEL_MAX;
    }
  else
    /* we only need to consider ideals of index >= j0, assuming the weight of
       an ideal cannot decrease (except when decreasing to zero when merged) */
    j0 = jmin[mat->cwmax];

  compute_weights_backend(mat, j0);

  if (jmin[0] == 0) /* jmin was not initialized */
    compute_jmin (mat, jmin);

  cpu = seconds () - cpu;
  wct = wct_seconds () - wct;
  print_timings ("   compute_weights took", cpu, wct);
  cpu_t[COMPUTE_W] += cpu;
  wct_t[COMPUTE_W] += wct;
}

/* computes the transposed matrix for columns of weight <= cwmax
 * (we only consider columns >= j0).
 *
 * Careful: "R" is overloaded.  This "R" is not the same as that of the main() function
 */
static void
compute_R (filter_matrix_t *mat, index_t j0)
{
  double cpu = seconds (), wct = wct_seconds ();

  index_t *Rp = mat->Rp;
  index_t *Rq = mat->Rq;
  index_t *Rqinv = mat->Rqinv;
  uint64_t nrows = mat->nrows;
  uint64_t ncols = mat->ncols;
  col_weight_t cwmax = mat->cwmax;

  /* compute the number of rows, the indices of the rows and the row pointers */

  int T = omp_get_max_threads();
  index_t tRnz[T];
  index_t tRn[T];
#pragma omp parallel
  {
          int T = omp_get_num_threads();
          int tid = omp_get_thread_num();
          index_t Rnz = 0;
          index_t Rn = 0;
          #pragma omp for schedule(static) nowait
          for (index_t j = j0; j < ncols; j++) {
              col_weight_t w = mat->wt[j];
              if (0 < w && w <= cwmax) {
                  Rnz += w;
                  Rn++;
              }
          }
          tRnz[tid] = Rnz;
          tRn[tid] = Rn;

#pragma omp barrier

          /* prefix-sum over the T threads (sequentially) */
#pragma omp single
          {
              index_t r = 0;
              index_t s = 0;
              for (int t = 0; t < T; t++) {
                  index_t w = tRnz[t];
                  index_t n = tRn[t];
                  tRnz[t] = r;
                  tRn[t] = s;
                  r += w;
                  s += n;
              }
              mat->Rn = s;
              Rp[s] = r; /* set the last row pointer */
          }
          Rnz = tRnz[tid];
          Rn = tRn[tid];

          #pragma omp for schedule(static) /* static is mandatory here */
          for (index_t j = j0; j < ncols; j++) {
              col_weight_t w = mat->wt[j];
              if (0 < w && w <= cwmax) {
                  Rq[j] = Rn;
                  Rqinv[Rn] = j;
                  Rnz += w;
                  Rp[Rn] = Rnz;
                  Rn++;
              }
          }
  } /* end parallel section */

  index_t Rn = mat->Rn;
  index_t Rnz = Rp[Rn];

  /* allocate variable-sized output (Rp is preallocated) */
  index_t *Ri = malloc_aligned (Rnz * sizeof(index_t), 64);
  mat->Ri = Ri;

  MAYBE_UNUSED double before_extraction = wct_seconds();

  /* dispatch entries */
  #pragma omp parallel for schedule(guided)
  for (index_t i = 0; i < nrows; i++) {
          if (mat->rows[i] == NULL)
                  continue; /* row was discarded */
          for (int k = matLengthRow(mat, i); k >= 1; k--) {
                  index_t j = matCell(mat, i, k);
                  if (j < j0)
                          break;
                  if (mat->wt[j] > cwmax)
                          continue;
                  index_t row = i;
                  index_t col = Rq[j];
                  uint64_t ptr;
                  #pragma omp atomic capture
                  ptr = --Rp[col];
                  Ri[ptr] = row;
          }
  }
  MAYBE_UNUSED double before_compression = wct_seconds();
  MAYBE_UNUSED double end_time = before_compression;

  cpu = seconds () - cpu;
  wct = wct_seconds () - wct;
  print_timings ("   compute_R took", cpu, wct);
  cpu_t[COMPUTE_R] += cpu;
  wct_t[COMPUTE_R] += wct;
}


static inline void
decrease_weight (filter_matrix_t *mat, index_t j)
{
  /* only decrease the weight if <= MERGE_LEVEL_MAX,
     since we saturate to MERGE_LEVEL_MAX+1 */
  if (mat->wt[j] <= MERGE_LEVEL_MAX) {
    /* update is enough, we do not need capture since we are not interested
       by the value of wt[j] */
    #pragma omp atomic update
    mat->wt[j]--;
#ifdef BIG_BROTHER_EXPENSIVE
    touched_columns[j] = 1;
#endif
  }
}

static inline void
increase_weight (filter_matrix_t *mat, index_t j)
{
  /* only increase the weight if <= MERGE_LEVEL_MAX,
     since we saturate to MERGE_LEVEL_MAX+1 */
  if (mat->wt[j] <= MERGE_LEVEL_MAX) {
    #pragma omp atomic update
    mat->wt[j]++;
#ifdef BIG_BROTHER_EXPENSIVE
    touched_columns[j] = 1;
#endif
  }
}

/* add row i2 to i1 (in place), and return the fill-in */
#ifndef FOR_DL

/* special code for factorization */
static int32_t
add_row (filter_matrix_t *L, index_t i1, index_t i2, MAYBE_UNUSED int32_t e1, MAYBE_UNUSED int32_t e2, buffer_struct_t *buf)
{
	typerow_t *r1 = L->rows[i1];
	typerow_t *r2 = L->rows[i2];
	index_t k1 = matLengthRow(L, i1);
	index_t k2 = matLengthRow(L, i2);
	index_t t1 = 1, t2 = 1;
	index_t t = 0;
        
        if (buf != NULL) {
                char s[MERGE_CHAR_MAX];
                int n = snprintf(s, MERGE_CHAR_MAX, "+%" PRId64 " %" PRId64 "\n", (uint64_t) i1, (uint64_t) i2);
                ASSERT(n < MERGE_CHAR_MAX);
                buffer_add(buf, s);
        }

#ifdef CANCEL
	#pragma omp atomic update
	cancel_rows ++;
#endif

	/* fast-track : don't precompute the size */
	typerow_t *sum = heap_alloc_row(L->heap, i1, k1 + k2);

	while (t1 <= k1 && t2 <= k2) {
		index_t j1 = r1[t1]; 
		index_t j2 = r2[t2];
		if (j1 == j2) {
			decrease_weight(L, j1);
			t1++;
			t2++;
		} else if (j1 < j2) {
			t++;
			sum[t] = j1;
			t1++;
		} else {
			increase_weight(L, j2);
			t++;
			sum[t] = j2;
			t2++;
		}
	}
	while (t1 <= k1)
	      sum[++t] = r1[t1++];
	while (t2 <= k2) {
	    increase_weight(L, r2[t2]);
	    sum[++t] = r2[t2++];
	}
	ASSERT(t <= k1 + k2);

#ifdef CANCEL
	int cancel = (t1 - 1) + (t2 - 1) - (t - 1);
	ASSERT_ALWAYS(cancel < CANCEL_MAX);
	#pragma omp atomic update
	cancel_cols[cancel] ++;
#endif

	heap_resize_last_row(L->heap, sum, t);
	heap_destroy_row(L->heap, L->rows[i1]);
	L->rows[i1] = sum;
	return t - k1; /* weight increase when replacing i1 by i1+i2 */
}

/* return the weight of the relation obtained when adding relations i1 and i2.
 * This is quite similar to add_row
 */
int
weightSum(const filter_matrix_t *L, index_t i1, index_t i2, MAYBE_UNUSED int32_t e1, MAYBE_UNUSED int32_t e2)
{
	typerow_t *r1 = L->rows[i1];
	typerow_t *r2 = L->rows[i2];
	index_t k1 = matLengthRow(L, i1);
	index_t k2 = matLengthRow(L, i2);
	index_t t1 = 1, t2 = 1;
	index_t t = 0;
	while (t1 <= k1 && t2 <= k2) {
		index_t j1 = r1[t1]; 
		index_t j2 = r2[t2];
		if (j1 == j2) {
			t1++;
			t2++;
		} else if (j1 < j2) {
			t++;
			t1++;
		} else {
			t++;
			t2++;
		}
	}
	t += (k1 - t1) + (k2 - t2);
	ASSERT(t <= k1 + k2);
	return t;
}


#else /* FOR_DL: j is the ideal to be merged */
#define INT32_MIN_64 (int64_t) INT32_MIN
#define INT32_MAX_64 (int64_t) INT32_MAX

/* we will multiply row i1 by e2, and row i2 by e1 */
static int32_t
add_row (filter_matrix_t *L, index_t i1, index_t i2,  int32_t e1, MAYBE_UNUSED int32_t e2, buffer_struct_t *buf)
{
	/* we always check e1 and e2 are not zero, in order to prevent from zero
     	 * exponents that would come from exponent overflows in previous merges 
     	 */
	ASSERT_ALWAYS (e1 != 0 && e2 != 0);
#ifdef CANCEL
	#pragma omp atomic update
	cancel_rows ++;
#endif

        if (buf != NULL) {
                char s[MERGE_CHAR_MAX];
                int n = snprintf(s, MERGE_CHAR_MAX, "*%" PRId64 " %" PRId32 " %" PRId64 " %" PRId32 "\n", (uint64_t) i1, e2, (uint64_t) i2, e1);
                ASSERT(n < MERGE_CHAR_MAX);
                buffer_add(buf, s);
        }

	/* now perform the real merge on L */
	index_t t1 = 1, t2 = 1, t = 0;
	typerow_t *r1 = L->rows[i1];
	typerow_t *r2 = L->rows[i2];
	index_t k1 = matLengthRow(L, i1);
	index_t k2 = matLengthRow(L, i2);
	typerow_t *sum = heap_alloc_row(L->heap, i1, k1 + k2);

	while (t1 <= k1 && t2 <= k2) {
		if (r1[t1].id == r2[t2].id) {
			/* as above, the exponent e below cannot overflow */
			int64_t e = (int64_t) e2 * (int64_t) r1[t1].e + (int64_t) e1 * (int64_t) r2[t2].e;
			if (e != 0) { /* exponents do not cancel */
				ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
				t++;
				setCell(sum, t, r1[t1].id, e);
			} else {
				decrease_weight (L, r1[t1].id);
			}
	  		t1++;
	  		t2++;
		} else if (r1[t1].id < r2[t2].id) {
			int64_t e = (int64_t) e2 * (int64_t) r1[t1].e;
			ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
			t++;
			setCell(sum, t, r1[t1].id, e);
			t1 ++;
		} else {
	  		int64_t e = (int64_t) e1 * (int64_t) r2[t2].e;
	  		ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
	  		t++;
	  		setCell(sum, t, r2[t2].id, e);
	  		increase_weight (L, r2[t2].id);
	  		t2 ++;
		}
	}
	while (t1 <= k1) {
		int64_t e = (int64_t) e2 * (int64_t) r1[t1].e;
		ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
		t++;
		setCell(sum, t, r1[t1].id, e);
		t1 ++;
	}
	while (t2 <= k2) {
		int64_t e = (int64_t) e1 * (int64_t) r2[t2].e;
		ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
		t++;
		setCell(sum, t, r2[t2].id, e);
		increase_weight (L, r2[t2].id);
		t2 ++;
	}
	ASSERT(t <= k1 + k2);

#ifdef CANCEL
	int cancel = (t1 - 1) + (t2 - 1) - (t - 1);
	ASSERT_ALWAYS(cancel < CANCEL_MAX);
	#pragma omp atomic update
	cancel_cols[cancel] ++;
#endif

	heap_resize_last_row(L->heap, sum, t);
	heap_destroy_row(L->heap, L->rows[i1]);
	L->rows[i1] = sum;
	return t - k1;
}

/* return the weight of the relation obtained when adding relations e1*i1 and e2*i2.
 * This is quite similar to add_row
 */
int
weightSum(const filter_matrix_t *L, index_t i1, index_t i2, MAYBE_UNUSED int32_t e1, MAYBE_UNUSED int32_t e2)
{
	ASSERT (e1 != 0 && e2 != 0);
	typerow_t *r1 = L->rows[i1];
	typerow_t *r2 = L->rows[i2];
	index_t k1 = matLengthRow(L, i1);
	index_t k2 = matLengthRow(L, i2);
	index_t t1 = 1, t2 = 1;
	index_t t = 0;

	while (t1 <= k1 && t2 <= k2) {
		if (r1[t1].id == r2[t2].id) {
			/* as above, the exponent e below cannot overflow */
			int64_t e = (int64_t) e2 * (int64_t) r1[t1].e + (int64_t) e1 * (int64_t) r2[t2].e;
			if (e != 0) 
				t++;   /* exponents do not cancel */
			t1++;
			t2++;
		} else if (r1[t1].id < r2[t2].id) {
			t++;
			t1++;
		} else {
			t++;
			t2++;
		}
	}
	t += (k1 - t1) + (k2 - t2);
	ASSERT(t <= k1 + k2);
	return t;
}

#endif

/* L is the "left" matrix (rows are relations-sets, columns are relations),
   and mat is the "merged" matrix (rows are relations-sets, columns are
   ideals). Return the (negative) fill-in.
*/
static int32_t
remove_row (filter_matrix_t *L, filter_matrix_t *mat, index_t i, buffer_struct_t *buf)
{
        if (buf != NULL) {
                char s[MERGE_CHAR_MAX];
                int n = snprintf(s, MERGE_CHAR_MAX, "-%" PRId64 "\n", (uint64_t) i);
                ASSERT(n < MERGE_CHAR_MAX);
                buffer_add(buf, s);
        }
        int32_t w = matLengthRow(mat, i);
        for (int k = 1; k <= w; k++)
                decrease_weight(mat, rowCell(mat->rows[i], k));
        heap_destroy_row(mat->heap, mat->rows[i]);
        mat->rows[i] = NULL;
        // replicate on L
        heap_destroy_row(L->heap, L->rows[i]);
        L->rows[i] = NULL;
        return -w;
}

#ifdef DEBUG
static void MAYBE_UNUSED
printRow (filter_matrix_t *mat, index_t i)
{
  if (mat->rows[i] == NULL)
    {
      printf ("row %u has been discarded\n", i);
      return;
    }
  int32_t k = matLengthRow (mat, i);
  printf ("%lu [%d]:", (unsigned long) i, k);
  for (int j = 1; j <= k; j++)
#ifndef FOR_DL
    printf (" %lu", (unsigned long) mat->rows[i][j]);
#else
    printf (" %lu^%d", (unsigned long) mat->rows[i][j].id, mat->rows[i][j].e);
#endif
  printf ("\n");
}
#endif

#if MERGE_STRATEGY == 1
#define BIAS 3
#else
#define BIAS (int32_t) ((col_weight_t) (-1))
#endif

/* classical cost: merge the row of smaller weight with the other ones,
   and return the merge cost (partially taking account of cancellations).
   Return:
   * 0 for 1-merges and 2-merges
   * a value >= 3 for all k-merges with k >= 3
   * here L might be mat (MERGE_STRATEGY=1) or L (MERGE_STRATEGY=2)
   */
static int32_t
merge_cost (filter_matrix_t *L, MAYBE_UNUSED filter_matrix_t *R,
            filter_matrix_t *mat, index_t j)
{
  index_t lo = mat->Rp[j];
  index_t hi = mat->Rp[j + 1];
  int32_t w = hi - lo; /* weight of the considered ideal */

  if (w <= 2)
    return 0; /* ensure all 1-and 2-merges are processed first
                 (they are done on the R matrix) */

  if (w > mat->cwmax) /* discard too heavy ideals for now */
    return INT32_MAX;

  /* find lightest row in L */
  index_t i = mat->Ri[lo];
  int32_t c, cmin = matLengthRow (L, i);
  for (index_t k = lo + 1; k < hi; k++)
    {
      i = mat->Ri[k];
      c = matLengthRow(L, i);
      if (c < cmin)
	cmin = c;
    }

  /* We return the original Markowitz cost, i.e., we optimize the weight
     increase of the matrix product M=L*R. */
#if MERGE_STRATEGY == 1
  ASSERT_ALWAYS(cmin >= 2);
  return (w - 1) * (cmin - 2) + BIAS;
#else
  return (w - 2) * cmin + BIAS - R->wt[j];
#endif
}


/* Rows i1, i2 in mat have an entry on column j; Set e1, e2 to the exponent */
static void locate_columns(MAYBE_UNUSED const filter_matrix_t *mat, 
	MAYBE_UNUSED index_t i1, MAYBE_UNUSED index_t i2, 
	MAYBE_UNUSED index_t j, MAYBE_UNUSED int32_t *e1, MAYBE_UNUSED int32_t *e2)
{
	#ifdef FOR_DL
	uint32_t k1 = matLengthRow (mat, i1);
	uint32_t k2 = matLengthRow (mat, i2);
	typerow_t *r1 = mat->rows[i1];
	typerow_t *r2 = mat->rows[i2];
	*e1 = 0;
	*e2 = 0;
		
	/* search by decreasing ideals as the ideal to be merged is likely large */
	for (int l = k1; l >= 1; l--)
		if (r1[l].id == j) {
			*e1 = r1[l].e;
			break;
		}
	for (int l = k2; l >= 1; l--)
		if (r2[l].id == j) {
			*e2 = r2[l].e;
			break;
		}
	int d = (int) gcd_int64 ((int64_t) *e1, (int64_t) *e2);
	*e1 /= -d;
	*e2 /= d;
	#endif
}

/* Perform the row additions given by the minimal spanning tree (stored in
   history[][]). Add in c the fill-in in mat. */
static int
addFatherToSons (index_t history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		 filter_matrix_t *L, filter_matrix_t *mat, int m, index_t *ind,
                 MAYBE_UNUSED index_t j, int *father, int *sons, int32_t *c, buffer_struct_t *buf)
{
	for (int i = m - 2; i >= 0; i--) {
		int s = father[i];
		int t = sons[i];
		if (i == 0) {
			history[i][1] = ind[s];
			ASSERT(s == 0);
		} else {
			history[i][1] = -(ind[s] + 1);
		}

		index_t i1 = ind[t];
		index_t i2 = ind[s];
		
		/* perform the operation on mat, and possibly "record" it in L 
		 * First determine the coefficients of the linear combination
		 */
		int32_t e1, e2;
		locate_columns(mat, i1, i2, j, &e1, &e2);

		*c += add_row (mat, i1, i2, e1, e2, buf);
		if (!twomerge_mode)
      			add_row (L, i1, i2, e1, e2, NULL);
		history[i][2] = i1;
		history[i][0] = 2;
	}
	return m - 2;
}


/* given an ideal of weight m, fills the m x m matrix A so that
   A[i][j] is the weight of the sum of the i-th and j-th rows
   containing the ideal, for 0 <= i, j < m */
void
fillRowAddMatrix(int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX],
		 const filter_matrix_t *mat, filter_matrix_t *L,
                 int m, index_t *ind, index_t j)
{
	/* A[i][i] is not used, thus we don't initialize it. */
	for (int u = 0; u < m; u++)
		for (int v = u+1; v < m; v++) {
			index_t i1 = ind[u];
			index_t i2 = ind[v];
			int32_t e1 = 0, e2 = 0;
			locate_columns(mat, i1, i2, j, &e1, &e2);
			A[u][v] = weightSum(L, i1, i2, e1, e2);
			A[v][u] = A[u][v];
		}
}


/* perform the merge described by the id-th row of R,
   computing the full spanning tree, and return the fill-in */
static int32_t
merge_do (filter_matrix_t *L, filter_matrix_t *R, filter_matrix_t *mat,
          index_t id, buffer_struct_t *buf)
{
	int32_t c = 0;
	index_t j = mat->Rqinv[id];
	index_t t = mat->Rp[id];
	int w = mat->Rp[id + 1] - t;

	ASSERT (1 <= w && w <= mat->cwmax);

        /* history: this eliminates column j */
	char s[MERGE_CHAR_MAX];
        int n = snprintf(s, MERGE_CHAR_MAX, "|%" PRId64 "\n", (uint64_t) mat->p[j]);
        ASSERT(n < MERGE_CHAR_MAX);
        buffer_add(buf, s);

	if (w == 1) {             // eliminate a singleton column
		index_signed_t i = mat->Ri[t]; /* only row containing j */        
		c = remove_row(L, mat, i, buf);
	} else {
		/* perform the real merge and output to history file */
		index_t *ind = mat->Ri + t;
		int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX];
#if MERGE_STRATEGY == 1
		fillRowAddMatrix (A, mat, mat, w, ind, j);           // target mat
#else
		fillRowAddMatrix (A, mat, L, w, ind, j);             // target L
#endif
		int start[MERGE_LEVEL_MAX], end[MERGE_LEVEL_MAX];
		minimalSpanningTree(start, end, w, A);
		index_t history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1];
		int _ MAYBE_UNUSED = addFatherToSons(history, L, mat, w, ind, j, start, end, &c, buf);
		c += remove_row(L, mat, ind[0], buf);
	}

	/* Update the total weight of R. We check the merged column has not a
	   saturated weight (which should not happen in practice). */
	ASSERT_ALWAYS (R->wt[j] < (col_weight_t) -1);
	if (twomerge_mode) {
		// we apply the merge directly to R
		#pragma omp atomic update
    		R->tot_weight += c;
	} else {
	        // we remove column j from R
    		/* Warning: R->wt[j] (the weight of column j after the 2-merges)
    		 * might differ from w (the current weight of column j). 
    		 */
		#pragma omp atomic update
		R->tot_weight -= R->wt[j];
	}
	return c;
}

/* accumulate in possible_merges all merges of (biased) cost <= cbound.
   possible_mergesL must be preallocated.
   possible_mergesL is a linear array and the merges appear by increasing cost.
   Returns the size of possible_merges. */
static int
compute_merges (index_t *possible_merges, MAYBE_UNUSED filter_matrix_t *L,
                MAYBE_UNUSED filter_matrix_t *R,
                filter_matrix_t *mat, int cbound)
{
  double cpu = seconds(), wct = wct_seconds();
  index_t Rn = mat->Rn;
  int * cost = malloc(Rn * sizeof(*cost));
  ASSERT_ALWAYS(cost != NULL);
  // int Lp[cbound + 2];  cost pointers

  /* compute the cost of all candidate merges */
  /* A dynamic schedule is needed here, since the columns of larger index have
     smaller weight, thus the load would not be evenly distributed with a
     static schedule. The value 128 was determined optimal experimentally
     on the RSA-512 benchmark with 32 threads, and is better than
     schedule(guided) for RSA-240 with 112 threads. */
  #pragma omp parallel for schedule(dynamic,128)
  for (index_t i = 0; i < Rn; i++)
#if MERGE_STRATEGY == 1
    cost[i] = merge_cost (mat, R, mat, i);
#else
    cost[i] = merge_cost (L, R, mat, i);
#endif

  int s;

  
  /* need an array that is visible to all threads in order to do the
   * prefix sum. Wish I knew another way.
   */
  int T = omp_get_max_threads();
  index_t count[T][cbound + 1];
  /* initialize array to zero */
#pragma omp for schedule(static)
  for (int t = 0; t < T; t++)
    memset(count+t, 0, (cbound + 1) * sizeof(index_t));

  /* Yet Another Bucket Sort (sigh): sort the candidate merges by cost. Check if worth parallelizing */
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    index_t *tcount = &count[tid][0];

#pragma omp for schedule(static)  // static is mandatory
    for (index_t i = 0; i < Rn; i++) {
      int c = cost[i];
      if (c <= cbound)
	tcount[c]++;
    }

         
#pragma omp single
    {
      /* prefix-sum */
      s = 0;
      for (int c = 0; c <= cbound; c++) {
	// Lp[c] = s;                     /* global row pointer in L */
	for (int t = 0; t < omp_get_num_threads(); t++) {
	  index_t w = count[t][c];       /* per-thread row pointer in L */
	  count[t][c] = s;
	  s += w;
	}
      }
    }

#pragma omp for schedule(static) // static is mandatory
    for (index_t i = 0; i < Rn; i++) {
      int c = cost[i];
      if (c > cbound)
	continue;
      possible_merges[tcount[c]++] = i;
    }
  } /* end parallel section */

  free(cost);

  double end = wct_seconds();
  double cpu2 = seconds() - cpu;
  double wct2 = end - wct;
  print_timings ("   compute_merges took", cpu2, wct2);
  cpu_t[COMPUTE_M] += cpu2;
  wct_t[COMPUTE_M] += wct2;
  return s;
}

/* return the number of merges applied */
static unsigned long
apply_merges (index_t *possible_merges, index_t total_merges,
              filter_matrix_t *L,  filter_matrix_t *R,
              filter_matrix_t *mat, buffer_struct_t *Buf)
{
  double cpu3 = seconds (), wct3 = wct_seconds ();
  char * busy_rows = malloc(mat->nrows * sizeof (char));
  memset (busy_rows, 0, mat->nrows * sizeof (char));

  unsigned long nmerges = 0;
  int64_t fill_in = 0;

  #pragma omp parallel reduction(+: fill_in, nmerges)
  {
    #pragma omp for schedule(guided)
    for (index_t it = 0; it < total_merges; it++) {
      index_t id = possible_merges[it];
      index_t lo = mat->Rp[id];
      index_t hi = mat->Rp[id + 1];
      int tid = omp_get_thread_num ();

      /* merge is possible if all its rows are "available" */
      int ok = 1;
      for (index_t k = lo; k < hi; k++) {
	index_t i = mat->Ri[k];
	if (busy_rows[i]) {
	  ok = 0;
	  break;
	}
      }
      if (ok) {
	  /* check again, since another thread might have reserved a row */
	  for (index_t k = lo; k < hi; k++) {
	    index_t i = mat->Ri[k];
	    char not_ok = 0;
	    /* we could use __sync_bool_compare_and_swap here,
	       but this is more portable and as efficient */
            #if defined(HAVE_OPENMP) && _OPENMP > 201107
	    /* the form of atomic capture below does not seem to be
	       recognized by OpenMP 3.5 (_OPENMP = 201107), see
	       https://cado-nfs-ci.loria.fr/ci/job/future-parallel-merge/job/compile-centos-6-i386/165 */
	    #pragma omp atomic capture
	    #else
	    #pragma omp critical
	    #endif
	    { not_ok = busy_rows[i]; busy_rows[i] = 1; }
	    if (not_ok)
	      {
		ok = 0;
		break;
	      }
	  }
      }
      if (ok) {
        fill_in += (uint64_t) merge_do (L, R, mat, id, &Buf[tid]);
        nmerges ++;
        ASSERT(hi - lo <= MERGE_LEVEL_MAX);
      }
    }  /* for */
  } /* parallel section */

  if (nmerges == 0 && total_merges > 0)
  {
    /* This can happen when we get a circular dependency between the row
       locks. In that case we simply apply the first merge.
       See https://gitlab.inria.fr/cado-nfs/cado-nfs/-/issues/30069. */
    index_t id = possible_merges[0];
    fill_in += merge_do(L, R, mat, id, Buf);
    nmerges ++;
    ASSERT(mat->Rp[id + 1] - mat->Rp[id] <= MERGE_LEVEL_MAX);
  }

  mat->tot_weight += fill_in;
  /* each merge decreases the number of rows and columns by one */
  mat->rem_nrows -= nmerges;
  mat->rem_ncols -= nmerges;

  double end = wct_seconds();
  free(busy_rows);
  cpu3 = seconds () - cpu3;
  wct3 = end - wct3;
  print_timings ("   apply_merges took", cpu3, wct3);
  cpu_t[APPLY_M] += cpu3;
  wct_t[APPLY_M] += wct3;
  return nmerges;
}

static double
average_density (filter_matrix_t *mat)
{
  return (double) mat->tot_weight / (double) mat->rem_nrows;
}

#ifdef DEBUG
/* duplicate the matrix, where the lines of mat_copy are
   not allocated individually by malloc, but all at once */
static void MAYBE_UNUSED
copy_matrix (filter_matrix_t *mat)
{
  unsigned long weight = mat->tot_weight;
  unsigned long s = weight + mat->nrows;
  index_t *T = malloc (s * sizeof (index_t));
  index_t *p = T;
  double cpu = seconds (), wct = wct_seconds ();
  for (index_t i = 0; i < mat->nrows; i++)
    {
      if (mat->rows[i] == NULL)
	{
	  p[0] = 0;
	  p ++;
	}
      else
	{
	  memcpy (p, mat->rows[i], (mat->rows[i][0] + 1) * sizeof (index_t));
	  p += mat->rows[i][0] + 1;
	}
    }
  print_timings ("   copy_matrix took", seconds () - cpu,
		 wct_seconds () - wct);
  ASSERT_ALWAYS(p == T + s);
  free (T);
}
#endif

#if 0
/* This function outputs the matrix in file 'out' in Sage format:
   M = matrix(...). Then to obtain a figure in Sage:
   sage: %runfile out.sage
   sage: M2 = matrix(RDF,512)
   sage: for i in range(512):
            for j in range(512):
               M2[i,j] = float(log(1+M[i,j]))
   sage: matrix_plot(M2,cmap='Greys')
   sage: matrix_plot(M2,cmap='Greys').save("mat.png")
*/
static void
output_matrix (filter_matrix_t *mat, char *out)
{
#define GREY_SIZE 512
    unsigned long grey[GREY_SIZE][GREY_SIZE];
    unsigned long hi = mat->nrows / GREY_SIZE;
    unsigned long hj = mat->ncols / GREY_SIZE;
    for (int i = 0; i < GREY_SIZE; i++)
      for (int j = 0; j < GREY_SIZE; j++)
	grey[i][j] = 0;
    for (index_t i = 0; i < mat->nrows; i++)
      {
	ASSERT_ALWAYS (mat->rows[i] != NULL);
	index_t ii = i / hi;
	if (ii >= GREY_SIZE)
	  ii = GREY_SIZE - 1;
	for (unsigned int k = 1; k <= matLengthRow(mat, i); k++)
	  {
	    index_t j = matCell(mat, i, k);
	    index_t jj = j / hj;
	    if (jj >= GREY_SIZE)
	      jj = GREY_SIZE - 1;
	    grey[ii][jj] ++;
	  }
      }
    FILE *fp = fopen (out, "w");
    fprintf (fp, "M=matrix([");
    for (int i = 0; i < GREY_SIZE; i++)
      {
	fprintf (fp, "[");
	int k = i;
	for (int j = 0; j < GREY_SIZE; j++)
	  {
	    fprintf (fp, "%lu", grey[k][j]);
	    if (j + 1 < GREY_SIZE)
	      fprintf (fp, ",");
	  }
	fprintf (fp, "]");
	if (i + 1 < GREY_SIZE)
	  fprintf (fp, ",");
      }
    fprintf (fp, "])\n");
    fclose (fp);
}
#endif

/*
 * This makes early verifications that our matrix isn't too large
 * compared to the types of the indices. Our goal is to bail out as
 * soon as we can.
 */
void sanity_check_matrix_sizes(filter_matrix_t * mat MAYBE_UNUSED)
{
#if (SIZEOF_INDEX == 4)
    if (mat->nrows >> 32)
    {
        fprintf (stderr, "Error, nrows = %" PRIu64 " larger than 2^32, please recompile with -DSIZEOF_INDEX=8\n", mat->nrows);
        exit (EXIT_FAILURE);
    }
    if (mat->ncols >> 32)
    {
        fprintf (stderr, "Error, ncols = %" PRIu64 " larger than 2^32, please recompile with -DSIZEOF_INDEX=8\n", mat->ncols);
        exit (EXIT_FAILURE);
    }
#endif
}

int
main (int argc, char *argv[])
{
    char *argv0 = argv[0];

    filter_matrix_t mat[1], L[1], R[1];
    FILE * history;

    int nthreads = 1, cbound_incr;
    uint32_t skip = DEFAULT_MERGE_SKIP;
    double target_density = DEFAULT_MERGE_TARGET_DENSITY;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    double tt;
    double cpu0 = seconds ();
    double wct0 = wct_seconds ();
    param_list pl;
    param_list_init (pl);
    declare_usage(pl);
    argv++,argc--;

    param_list_configure_switch (pl, "-v", &merge_verbose);
    param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    if (argc == 0)
      usage (pl, argv0);

    for( ; argc ; ) {
      if (param_list_update_cmdline(pl, &argc, &argv)) continue;
      fprintf (stderr, "Unknown option: %s\n", argv[0]);
      usage (pl, argv0);
    }
    /* print command-line arguments */
    verbose_interpret_parameters (pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    const char *purgedname = param_list_lookup_string (pl, "mat");
    const char *outname = param_list_lookup_string (pl, "out");
    const char *path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");

    param_list_parse_int (pl, "t", &nthreads);
#ifdef HAVE_OPENMP
    omp_set_num_threads (nthreads);
#endif

    if (param_list_parse_int (pl, "incr", &cbound_incr) == 0)
      cbound_incr = CBOUND_INCR_DEFAULT;

    param_list_parse_uint (pl, "skip", &skip);

    param_list_parse_double (pl, "target_density", &target_density);

    /* Some checks on command line arguments */
    if (param_list_warn_unused(pl))
    {
      fprintf(stderr, "Error, unused parameters are given\n");
      usage(pl, argv0);
    }

    if (purgedname == NULL)
    {
      fprintf(stderr, "Error, missing -mat command line argument\n");
      usage (pl, argv0);
    }
    if (outname == NULL)
    {
      fprintf(stderr, "Error, missing -out command line argument\n");
      usage (pl, argv0);
    }

    set_antebuffer_path (argv0, path_antebuffer);

    history = fopen_maybe_compressed (outname, "w");
    ASSERT_ALWAYS(history != NULL);

    /* some explanation about the history file */
    fprintf (history, "# Please refer to filter/README.fileformats for the syntax of this file.\n");

    /* Read number of rows and cols on first line of purged file */
    purgedfile_read_firstline (purgedname, &(mat->nrows), &(mat->ncols));

    sanity_check_matrix_sizes(mat);

    /* initialize the matrix structure */
    initMat (mat, skip); // merge_replay_matrix.c

    /* Set L to the identity matrix. This matrix will contain the combinations
       of rows of the initial matrix M, so that the output matrix M' of merge
       equals to L*M (up to removed rows and columns).
    */
    L->nrows = L->ncols = mat->nrows;
    initMat (L, 0); // no skip in L, since skip only concerns columns
    for (index_t i = 0; i < L->nrows; i++)
    {
      L->rows[i] = heap_alloc_row(L->heap, i, 2);  // merge_heap.c
      setCell(L->rows[i], 0, 1, 1); // only one non-zero element in row
      setCell(L->rows[i], 1, i, 1); // L[i,i] = 1
    }

    /* Allocate the wt field of R (we don't use the other fields). */
    R->wt = malloc (mat->ncols * sizeof (col_weight_t));

    /* Read all rels and fill-in the mat structure */
    tt = seconds ();
    filter_matrix_read (mat, purgedname);
    printf ("Time for filter_matrix_read: %2.2lfs\n", seconds () - tt);

    check_matrix (mat);

    buffer_struct_t *Buf = buffer_init (nthreads);

    double cpu_after_read = seconds ();
    double wct_after_read = wct_seconds ();

    /* jmin[w] for 1 <= w <= MERGE_LEVEL_MAX is the smallest column of weight w
       at beginning. We set jmin[0] to 0 to tell that jmin[] was not
       initialized. */
    index_t jmin[MERGE_LEVEL_MAX + 1] = {0,};

    compute_weights (mat, jmin);
    /* Initially, R=mat, thus they have the same weights. */
    compute_weights_R (R, mat);
    R->tot_weight = mat->tot_weight;

    recompress (R, mat, jmin);

    // output_matrix (mat, "out.sage");

    /* Allocate the transposed matrix R in CSR format. Since Rp is of fixed
       size, we allocate it for once. However, the size of Ri will vary from
       step to step. */
    mat->Rp = malloc ((mat->ncols + 1) * sizeof (index_t));
    mat->Ri = NULL;
    mat->Rq = malloc (mat->ncols * sizeof (index_t));
    mat->Rqinv = malloc (mat->ncols * sizeof (index_t));

    printf ("Using MERGE_LEVEL_MAX=%d, cbound_incr=%d",
	    MERGE_LEVEL_MAX, cbound_incr);
    printf (", PAGE_DATA_SIZE=%d", heap_config_get_PAGE_DATA_SIZE());
#ifdef HAVE_OPENMP
    /* https://stackoverflow.com/questions/38281448/how-to-check-the-version-of-openmp-on-windows
       201511 is OpenMP 4.5 */
    printf (", OpenMP %d", _OPENMP);
#endif
    printf ("\n");

    printf ("N=%" PRIu64 " W=%" PRIu64 " W/N=%.2f cpu=%.1fs wct=%.1fs mem=%luM\n",
	    mat->rem_nrows, mat->tot_weight, average_density (mat),
	    seconds () - cpu0, wct_seconds () - wct0,
	    PeakMemusage () >> 10);

    fflush (stdout);

    mat->cwmax = 2;

#if defined(DEBUG) && defined(FOR_DL)
    /* compute the minimum/maximum coefficients */
    int32_t min_exp = 0, max_exp = 0;
    for (index_t i = 0; i < mat->nrows; i++)
      if (mat->rows[i] != NULL)
	{
	  ideal_merge_t *ri = mat->rows[i];
	  for (unsigned int k = 1; k <= matLengthRow (mat, i); k++)
	    {
	      int32_t e = ri[k].e;
	      if (e < min_exp)
		min_exp = e;
	      if (e > max_exp)
		max_exp = e;
	    }
	}
    printf ("min_exp=%d max_exp=%d\n", min_exp, max_exp);
#endif

    unsigned long lastN, lastW;
    double lastWoverN;
    int cbound = BIAS; /* bound for the (biased) cost of merges to apply */
    int merge_pass = 0;

    /****** begin main loop ******/
    while (1) {
	double cpu1 = seconds (), wct1 = wct_seconds ();
	merge_pass++;

    if (merge_pass == 2 || mat->cwmax > 2) {
            double cpu8 = seconds (), wct8 = wct_seconds ();
            heap_garbage_collection(mat->heap, mat->rows);
            cpu8 = seconds () - cpu8;
            wct8 = wct_seconds () - wct8;
            print_timings ("   GC took", cpu8, wct8);
            cpu_t[GC] += cpu8;
            wct_t[GC] += wct8;
    }

	/* Once cwmax >= 3, at each pass, we increase cbound to allow more
	   merges. If one decreases cbound_incr, the final matrix will be
	   smaller, but merge will take more time.
	   If one increases cbound_incr, merge will be faster, but the final
	   matrix will be larger. */
	if (mat->cwmax > 2)
		cbound += cbound_incr;

	lastN = mat->rem_nrows;
	lastW = mat->tot_weight;
	lastWoverN = (double) lastW / (double) lastN;

	#ifdef TRACE_J
	for (index_t i = 0; i < mat->ncols; i++) {
		if (mat->rows[i] == NULL)
			continue;
		for (index_t k = 1; k <= matLengthRow(mat, i); k++)
                  if (rowCell(mat->rows[i],k) == TRACE_J)
                    printf ("ideal %d in row %lu\n", TRACE_J, (unsigned long) i);
	}
	#endif

	compute_R (mat, jmin[mat->cwmax]);

	index_t *possible_merges = malloc(mat->Rn * sizeof(index_t));

	index_t n_possible_merges = compute_merges(possible_merges, L, R, mat,
                                                   cbound);

	unsigned long nmerges = apply_merges (possible_merges,
                                            n_possible_merges, L, R, mat, Buf);

	buffer_flush (Buf, nthreads, history);
	free(possible_merges);

	free_aligned (mat->Ri);

        if (nmerges == 0 && n_possible_merges > 0)
          {
            fprintf (stderr, "Error, no merge done while n_possible_merges > 0\n");
            fprintf (stderr, "Please check the entries in your purged file are sorted\n");
            exit (EXIT_FAILURE);
          }

	/* settings for next pass */
  	if (mat->cwmax == 2) { /* we first process all 2-merges */
		if (nmerges == 0) {
                        twomerge_mode = 0;    // stop special treatment of the initial run of 2-merges
                        fprintf(history, "! End of two-merges\n");
                        fflush(history);
                        compute_weights_R (R, mat);
			mat->cwmax++;
                }
	} else {
#define INCREASE_RATIO 0.005
          /* we only increase cwmax when nmerges is below a given proportion
             of the matrix size */
          if (mat->cwmax < MERGE_LEVEL_MAX &&
              (nmerges < INCREASE_RATIO * (double) mat->ncols))
              mat->cwmax ++;
	}

	if (mat->rem_ncols < 0.66 * mat->ncols) {
	  static int recompress_pass = 0;
	  printf("============== Recompress %d ==============\n", ++recompress_pass);
	  recompress (R, mat, jmin);
	}

	cpu1 = seconds () - cpu1;
	wct1 = wct_seconds () - wct1;
	print_timings ("   pass took", cpu1, wct1);
	cpu_t[PASS] += cpu1;
	wct_t[PASS] += wct1;

	/* estimate current average fill-in */
	double av_fill_in = ((double) mat->tot_weight - (double) lastW)
	  / (double) (lastN - mat->rem_nrows);

	printf ("N=%" PRIu64 " W=%" PRIu64 " (%.0fMB) W/N=%.2f fill-in=%.2f cpu=%.1fs wct=%.1fs mem=%luM [pass=%d,cwmax=%d]\n",
		mat->rem_nrows, mat->tot_weight,
		9.5367431640625e-07 * (mat->rem_nrows + mat->tot_weight) * sizeof(index_t),
		(double) mat->tot_weight / (double) mat->rem_nrows, av_fill_in,
		seconds () - cpu0, wct_seconds () - wct0,
		PeakMemusage () >> 10, merge_pass, mat->cwmax);
	fflush (stdout);

	if (average_density (mat) >= target_density)
          break;

        /* With small cbound_incr, in particular cbound_incr=1,
           we might have zero potential merge when cbound is small,
           thus we stop only when cbound > cwmax^2 (the cost of a
           merge being proportional to the square of the column weight). */
	if (nmerges == 0 && mat->cwmax == MERGE_LEVEL_MAX &&
            cbound > mat->cwmax * mat->cwmax)
          break;
    }
    /****** end main loop ******/

#if defined(DEBUG) && defined(FOR_DL)
    min_exp = 0; max_exp = 0;
    for (index_t i = 0; i < mat->nrows; i++)
      if (mat->rows[i] != NULL)
	{
	  for (unsigned int k = 1; k <= matLengthRow (mat, i); k++)
	    {
	      int32_t e = mat->rows[i][k].e;
	      if (e < min_exp)
		min_exp = e;
	      if (e > max_exp)
		max_exp = e;
	    }
	}
    printf ("min_exp=%d max_exp=%d\n", min_exp, max_exp);
#endif

    fclose_maybe_compressed (history, outname);

    if (average_density (mat) > target_density)
      {
	/* estimate N for W/N = target_density, assuming W/N = a*N + b */
	unsigned long N = mat->rem_nrows;
	double WoverN = (double) mat->tot_weight / (double) N;
	double a = (lastWoverN - WoverN) / (double) (lastN - N);
	double b = WoverN - a * (double) N;
	/* we want target_density = a*N_target + b */
	printf ("Estimated N=%" PRIu64 " for W/N=%.2f\n",
		(uint64_t) ((target_density - b) / a), target_density);
      }

    print_timings ("compute_weights:", cpu_t[COMPUTE_W], wct_t[COMPUTE_W]);
    print_timings ("compute_R      :", cpu_t[COMPUTE_R], wct_t[COMPUTE_R]);
    print_timings ("compute_merges :", cpu_t[COMPUTE_M], wct_t[COMPUTE_M]);
    print_timings ("apply_merges   :", cpu_t[APPLY_M], wct_t[APPLY_M]);
    print_timings ("pass           :", cpu_t[PASS], wct_t[PASS]);
    print_timings ("recompress     :", cpu_t[RECOMPRESS], wct_t[RECOMPRESS]);
    print_timings ("buffer_flush   :", cpu_t[FLUSH], wct_t[FLUSH]);
    print_timings ("garbage_coll   :", cpu_t[GC], wct_t[GC]);

    printf ("Final matrix has N=%" PRIu64 " nc=%" PRIu64 " (%" PRIu64
	    ") W=%" PRIu64 "\n", mat->rem_nrows, mat->rem_ncols,
	    mat->rem_nrows - mat->rem_ncols, mat->tot_weight);
    fflush (stdout);

    /* stats on L (mostly for debugging for now) */
    uint64_t Lweight = 0;
    uint64_t Ln = 0;
    for (index_t i = 0; i < L->nrows; i++)
        if (L->rows[i] != NULL) {
                Ln += 1;
                Lweight += rowLength(L->rows, i);
        }

    printf ("L has N=%" PRIu64 " W=%" PRIu64 "\n", Ln, Lweight);
    printf ("R has N=%" PRIu64 " W=%" PRIu64 "\n",
            mat->rem_ncols, R->tot_weight);
    printf ("Theoretical linalg cost N*(|L|+|R|)=%.2e\n",
            (double) Ln * ((double) Lweight + (double) R->tot_weight));
    fflush (stdout);

/*
    FILE *f = fopen("Ldump.txt", "w");
    for (index_t i = 0; i < L->nrows; i++) 
        if (L->rows[i] != NULL) {
                for (index_t j = 0; j <= rowLength(L->rows, i); j++)
                        fprintf(f, "%d ", rowCell(L->rows[i], j));
                fprintf(f, "\n");
        }
    fclose(f);
*/
    printf ("Before cleaning memory:\n");
    print_timing_and_memory (stdout, cpu_after_read, wct_after_read);

    buffer_clear (Buf, nthreads);

    free (mat->p);
    free (mat->Rp);
    free (mat->Rq);
    free (mat->Rqinv);

    clearMat (mat);
    clearMat (L);
    free (R->wt);

    param_list_clear (pl);

    printf ("After cleaning memory:\n");
    print_timing_and_memory (stdout, cpu_after_read, wct_after_read);

    /* print total time and memory (including reading the input matrix,
       initializing and free-ing all data) */
    print_timing_and_memory (stdout, cpu0, wct0);

#ifdef CANCEL
    unsigned long tot_cancel = 0;
    printf ("cancel_rows=%lu\n", cancel_rows);
    for (int i = 0; i < CANCEL_MAX; i++)
      if (cancel_cols[i] != 0)
	{
	  tot_cancel += cancel_cols[i] * i;
	  printf ("cancel_cols[%d]=%lu\n", i, cancel_cols[i]);
	}
    printf ("tot_cancel=%lu\n", tot_cancel);
#endif

    return 0;
}
