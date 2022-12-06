/* TODO:
 * modify merge to give a termination criteria in the double matrix case:
   we suggest to use the average density of the left matrix L
   (-target_density 100 would mean that we stop when the average density
   of the rows of L reaches 100)
 * maybe modify the criteria use in merge to choose the next merge to perform.
   Currently we use the merge that minimizes the fill-in in M = L*R, but in
   the double matrix case, we want to minimize the fill-in in L.
 */

/* replay --- replaying history of merges to build the sparse matrices in
              Kleinjung's "double matrix" idea.

Copyright 2008-2022 Francois Morain, Emmanuel Thome, Paul Zimmermann,
          Cyril Bouvier, Pierrick Gaudry, Charles Bouillaguet

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

#include "cado.h"		// IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MINGW
#include <fcntl.h>		/* for _O_BINARY */
#endif
#include <string.h>
#include <inttypes.h>		// for PRIu64, PRIu32, PRIx64
#include <stdint.h>		// for uint64_t, uint32_t, UINT32_MAX
#include "purgedfile.h"		// for purgedfile_read_firstline
#include "typedefs.h"		// for index_t, ideal_merge_t, index_signed_t
#include "filter_config.h"
#include "filter_io.h"		// earlyparsed_relation_ptr
#include "fix-endianness.h"	// fwrite32_little
#include "gzip.h"		// fopen_maybe_compressed
#include "misc.h"		// derived_filename, has_suffix
#include "params.h"		// param_list_parse_*
#include "sparse.h"
#include "stats.h"		// stats_data_t
#include "timing.h"		// seconds
#include "verbose.h"		// verbose_decl_usage
#include "portability.h"	// strdup  // IWYU pragma: keep
#include "macros.h"
#include "merge_heap.h"
#include "gcd.h"

#define DEBUG 0

/* purged matrix */
typerow_t **rows;
index_t nrows;          // #rows for the purged matrix
index_t ncols;          // #cols for the purged matrix
index_t nrows_small;    // #rows for the purged matrix once 2-merges are performed        

/* global state */
index_t *column_info;
char *scratch;


/* Save a sparse matrix to the filesystem. Write column and row weights in
   companion files. If skip > 0, also create a "dense" part in the same format. 
   rows_size == #entries in rows[].
   nrows     == #non-NULL entries in rows[]
*/
static unsigned long flushSparse(const char *sparsename, typerow_t ** rows,
                           index_t rows_size, index_t nrows, index_t ncols, index_t skip, int bin)
{
	double tt = seconds();
	printf("Sparse submatrix: nrows=%" PRIu64 " ncols=%" PRIu64 "\n",
		(uint64_t) nrows, (uint64_t) ncols);
	printf("Writing sparse representation to %s\n", sparsename);
	fflush(stdout);

#ifdef FOR_DL
	ASSERT_ALWAYS(skip == 0);
#endif
	const struct {
		const char *ext;
		const char *smat;
		const char *srw;
		const char *scw;
		const char *dmat;
		const char *drw;
		const char *dcw;
	} suffixes[2] = {
        {
          .ext = ".txt",
          .smat = "txt",
          .srw = "rw.txt",
          .scw = "cw.txt",
          .dmat = "dense.txt",
          .drw = "dense.rw.txt",
          .dcw = "dense.cw.txt",
        },
        {
          .ext = ".bin",
          .smat = "bin",
          .srw = "rw.bin",
          .scw = "cw.bin",
          .dmat = "dense.bin",
          .drw = "dense.rw.bin",
          .dcw = "dense.cw.bin",
        },
        }, * suf = &(suffixes[bin]);
	unsigned long W = 0;  /* sparse weight */
	unsigned long DW = 0; /* dense weight */
	char *zip = NULL;
	index_t *weights = malloc(ncols * sizeof(index_t));
	ASSERT_ALWAYS(weights != NULL);
	memset(weights, 0, ncols * sizeof(index_t));

	/* [...] the 'b' is ignored on all POSIX conforming systems, including Linux.
	 * (Other systems may treat text files and binary files
	 * differently, and adding the 'b' may  be a good idea if you do
	 * I/O to a binary file and expect that your program may be ported
	 * to non-UNIX environments.)
	 * ---- fopen man page */
	char wmode[3] = "wb";

	/* setup filenames */
	char *base = strdup(sparsename);
	if (has_suffix(base, suf->ext)) { /* strip suffix if given */
		base[strlen(base) - 4] = '\0';
	}
	char *smatname = NULL;
	FILE *smatfile = NULL;
	char *srwname = NULL;
	FILE *srwfile = NULL;
	char *scwname = NULL;
	FILE *scwfile = NULL;
	char *dmatname = NULL;
	FILE *dmatfile = NULL;
	char *drwname = NULL;
	FILE *drwfile = NULL;
	char *dcwname = NULL;
	FILE *dcwfile = NULL;

	smatname = derived_filename(base, suf->smat, zip);
	srwname = derived_filename(base, suf->srw, zip);
	scwname = derived_filename(base, suf->scw, zip);
	smatfile = fopen(smatname, wmode);
	srwfile = fopen(srwname, wmode);
        if (!bin) fprintf(smatfile, "%" PRIu64 " %" PRIu64 "\n",
                          (uint64_t) nrows,
                          (uint64_t) ncols - skip);

	if (skip) {
		/* arrange so that we don't get file names like .sparse.dense */
		char *dbase = strdup(base);
		char *tmp = strstr(dbase, ".sparse");
		if (tmp)
			memmove(tmp, tmp + 7, strlen(tmp + 7) + 1);
		dmatname = derived_filename(dbase, suf->dmat, zip);
		drwname = derived_filename(dbase, suf->drw, zip);
		dcwname = derived_filename(dbase, suf->dcw, zip);
		dmatfile = fopen(dmatname, wmode);
		drwfile = fopen(drwname, wmode);
                if (!bin)
                  fprintf(dmatfile, "%" PRIu64 " %" PRIu64 "\n",
                          (uint64_t) nrows, (uint64_t) skip);
		free(dbase);
	}

	for (index_t i = 0; i < rows_size; i++) {
		if (rows[i] == NULL) {
			continue; /* row has been deleted */
		} else {
			uint32_t dw = 0;
			uint32_t sw = 0;
			for (index_t j = 1; j <= rowLength(rows, i); j++) {
				if (rowCell(rows[i], j) < skip) {
					dw++;
					DW++;
				} else {
					sw++;
					W++;
				}
			}
                        if (bin)
                        {
                          fwrite32_little(&sw, 1, smatfile);
                          if (srwfile) fwrite32_little(&sw, 1, srwfile);
                          if (skip) fwrite32_little(&dw, 1, dmatfile);
                          if (skip) fwrite32_little(&dw, 1, drwfile);
                        } else {
                          fprintf(smatfile, "%" PRIu32 "", sw);
                          if (srwfile) fprintf(srwfile, "%" PRIu32 "\n", sw);
                          if (skip) fprintf(dmatfile, "%" PRIu32 "", dw);
                          if (skip) fprintf(drwfile, "%" PRIu32 "\n", dw);
                        }

			for (index_t j = 1; j <= rowLength(rows, i); j++) {
				ASSERT_ALWAYS(rowCell(rows[i], j) <= (index_t) UINT32_MAX);
				uint32_t x = rowCell(rows[i], j);
				if (srwfile)
					weights[x]++;
				if (x < skip) {
                                  ASSERT_ALWAYS(skip);
                                  if (bin)
                                    fwrite32_little(&x, 1, dmatfile);
                                  else
                                    fprintf(dmatfile, " %" PRIu32 "", x);
				} else {
					x -= skip;
                                        if (bin) {
                                          fwrite32_little(&x, 1, smatfile);
#ifdef FOR_DL
					/* exponents are always int32_t */
					uint32_t e = (uint32_t) rows[i][j].e;
					fwrite32_little(&e, 1, smatfile);
#endif
                                        } else {
                                          fprintf(smatfile, " %" PRIu32 "", x);
#ifdef FOR_DL
                                          fprintf(smatfile, ":%d", rows[i][j].e);
#endif
                                        }
				}
			}
		}
                if (!bin) {
                  fprintf(smatfile, "\n");
                  if (skip) fprintf(dmatfile, "\n");
                }
	}

	fclose(smatfile);
	if (srwfile)
		fclose(srwfile);

	if (skip) {
		printf("%lu coeffs (out of %lu total) put into %s (%.1f%%)\n",
			DW, DW + W, dmatname, 100.0 * (double)DW / (DW + W + (DW == 0 && W == 0)));
		fflush(stdout);
		fclose(dmatfile);
		fclose(drwfile);
		dcwfile = fopen(dcwname, wmode);
		for (index_t j = 0; j < skip; j++) {
			ASSERT_ALWAYS(weights[j] <= (index_t) UINT32_MAX);
			uint32_t x = weights[j];
			fwrite32_little(&x, 1, dcwfile);
		}
		fclose(dcwfile);
	}
	scwfile = fopen(scwname, wmode);
	for (index_t j = skip; j < ncols; j++) {
		ASSERT_ALWAYS(weights[j] <= (index_t) UINT32_MAX);
		uint32_t x = weights[j];
		fwrite32_little(&x, 1, scwfile);
	}
	fclose(scwfile);
	free(smatname);
	free(srwname);
	free(scwname);
	free(dmatname);
	free(drwname);
	free(dcwname);
	free(base);
	free(weights);

	printf("# Writing matrix took %.1lfs\n", seconds() - tt);
	printf("# Weight of the sparse submatrix: %lu\n", W);
	fflush(stdout);

	return W;
}

#define STRLENMAX 2048

// i1 += i2
// j is the index of the column that is used for
// pivoting in the case of DL. Then, the operation is
//   i1 = e2*i1 + e1*i2
// where e1 and e2 are adjusted so that the j-th column is zero in i1.

/******************************************************************************/


/*  Set eliminated_columns[j] == 1 if column j is eliminated */
static void
preread_history(const char *hisname)
{
	printf("Reading history file %s (1st pass)\n", hisname);
	fflush(stdout);

	FILE * hisfile = fopen_maybe_compressed(hisname, "r");
	ASSERT_ALWAYS(hisfile != NULL);

	/* allocate auxiliary data (per column) */
	column_info = malloc(sizeof(*column_info) * ncols);
	ASSERT_ALWAYS(column_info != NULL);
	for (index_t i = 0; i < ncols; i++)
		column_info[i] = 0;

	uint64_t addread = 0;
	char str[STRLENMAX];
	stats_data_t stats;	/* struct for printing progress */
	stats_init(stats, stdout, &addread, 23, "Read", "row additions", "", "lines");
	uint64_t n_elim = 0;
	while (fgets(str, STRLENMAX, hisfile)) {
		if (str[0] == '#')
			continue;
		addread++;
		if (stats_test_progress(stats))
			stats_print_progress(stats, addread, 0, 0, 0);
		if (str[strlen(str) - 1] != '\n') {
			fprintf(stderr, "Gasp: not a complete line!");
			fprintf(stderr, " I stop reading and go to the next phase\n");
			break;
		}
		index_t j;
		index_signed_t ind[MERGE_LEVEL_MAX];
		int MAYBE_UNUSED _ = parse_hisfile_line(ind, str, &j);   // in sparse.c, mutualized with "normal" replay
		if (column_info[j] != 1)
			n_elim += 1;
		column_info[j] = 1;         // column has been eliminated
	}
	stats_print_progress(stats, addread, 0, 0, 1);
	fclose_maybe_compressed(hisfile, hisname);
	printf("%" PRId64 " eliminated columns\n", n_elim);
}


/*
 * The following code is very similar to replay.c
 * It loads the whole matrix in memory.
 */

void * read_purged_row (void MAYBE_UNUSED *context_data, earlyparsed_relation_ptr rel)
{
	typerow_t buf[UMAX(weight_t)];

	// 1st pass, set scratch
	for (unsigned int j = 0; j < rel->nb; j++) {
		index_t h = rel->primes[j].h;
		if (column_info[h] == 1)          // column was eliminated
			continue;
		scratch[h] ^= 1;
	}

	// 2nd pass, read and reset scratch
	unsigned int nb = 0;
	for (unsigned int j = 0; j < rel->nb; j++) {
		index_t h = rel->primes[j].h;
		if (column_info[h] == 1)          // column was eliminated
			continue;
		if (scratch[h] == 0)
			continue; // ideal appeared an even number of times
		scratch[h] = 0;             // reset scratch
		column_info[h] = 2;         // column is not empty
		nb += 1;
		#ifdef FOR_DL
			exponent_t e = rel->primes[j].e;
			buf[nb] = (ideal_merge_t) {.id = h, .e = e};
		#else
			ASSERT_ALWAYS (rel->primes[j].e == 1);
			buf[nb] = h;
		#endif
	}
	#ifdef FOR_DL
		buf[0].id = nb;
	#else
		buf[0] = nb;
	#endif
	
	/* required because of add_rows (2-merges) on the relations */
	qsort (&(buf[1]), nb, sizeof(typerow_t), cmp_typerow_t);

	rows[rel->num] = heap_alloc_row(rel->num, nb);  // in merge_heap.c
	compressRow (rows[rel->num], buf, nb);              // in sparse.c
  	return NULL;
}

/* update rows, column_info. 
 * Set column_info[j] == 2 if column j is non-empty 
 * Set column_info[j] == 0 otherwise 
 */
static void
read_purgedfile (const char* purgedname)
{
	printf("Reading purged matrix from %s\n", purgedname);
	fflush(stdout);

	/* allocate purged matrix */
	rows = malloc(nrows * sizeof(*rows));
	ASSERT_ALWAYS(rows != NULL);
	for (index_t i = 0; i < nrows; i++)
		rows[i] = NULL;

	char *fic[2] = {(char *) purgedname, NULL};
	scratch = malloc(ncols * sizeof(*scratch));
	ASSERT_ALWAYS(scratch != NULL);
	for (index_t i = 0; i < ncols; i++)
		scratch[i] = 0;

	index_t nread = filter_rels(fic, (filter_rels_callback_t) &read_purged_row, 
				NULL, EARLYPARSE_NEED_INDEX, NULL, NULL);
	ASSERT_ALWAYS (nread == nrows);
	free(scratch);
	/* here: column_info[j] == 1   <====>   column has been eliminated
		 column_info[j] == 2   <====>   column is non-empty (not eliminated)
		 column_info[j] == 0   <====>   column is empty (not eliminated) */
}


/* COPIED-PASTED from merge.c. REFACTORING PLAN: take rows[] and weights[] as
   arugments, update weights only if not NULL */

#ifndef FOR_DL
/* special code for factorization */
static void
add_row (typerow_t **rows, index_t i1, index_t i2, MAYBE_UNUSED index_t j)
{
	typerow_t *r1 = rows[i1];
  	typerow_t *r2 = rows[i2];
	index_t k1 = rowLength(rows, i1);
	index_t k2 = rowLength(rows, i2);
	index_t t1 = 1;            // index in rows[i1]
	index_t t2 = 1;            // index in rows[i2]
	index_t t = 0;             // index in the sum

	/* fast-track : don't precompute the size */
	typerow_t *sum = heap_alloc_row(i1, k1 + k2);

	while (t1 <= k1 && t2 <= k2) {
		if (r1[t1] == r2[t2]) {
			t1 += 1;  // cancellation
			t2 += 1;
		} else if (r1[t1] < r2[t2]) {
			t += 1;
			sum[t] = r1[t1];
			t1 += 1;
		} else {
			t += 1;
			sum[t] = r2[t2];
			t2 += 1;
		}
	}
	while (t1 <= k1) {
		t += 1;
		sum[t] = r1[t1];
		t1 += 1;
	}
	while (t2 <= k2) {
		t += 1;
		sum[t] = r2[t2];
		t2 += 1;
	}
        /* In the double-matrix code, we don't necessarily have
           cancellations. */
	ASSERT(t <= k1 + k2);
	heap_resize_last_row(sum, t);
	heap_destroy_row(r1);
	rows[i1] = sum;
	return;
}
#else /* FOR_DL: j is the ideal to be merged */
#define INT32_MIN_64 (int64_t) INT32_MIN
#define INT32_MAX_64 (int64_t) INT32_MAX

static void
add_row (typerow_t **rows, index_t i1, index_t i2, index_t j)
{
  /* first look for the exponents of j in i1 and i2 */
 	typerow_t *r1 = rows[i1];
  	typerow_t *r2 = rows[i2];
 	index_t k1 = rowLength(rows, i1);
	index_t k2 = rowLength(rows, i2);
	index_t t1 = 1;            // index in rows[i1]
	index_t t2 = 1;            // index in rows[i2]
	index_t t = 0;             // index in the sum

  	int32_t e1 = 0;
  	int32_t e2 = 0;

  /* search by decreasing ideals as the ideal to be merged is likely large */
  for (int l = k1; l >= 1; l--)
    if (r1[l].id == j) {
	e1 = r1[l].e;
	break;
      }
  for (int l = k2; l >= 1; l--)
    if (r2[l].id == j) {
	e2 = r2[l].e;
	break;
      }

  /* we always check that e1 and e2 are not zero, in order to prevent from zero
     exponents that would come from exponent overflows in previous merges */
  ASSERT_ALWAYS (e1 != 0 && e2 != 0);

  int d = (int) gcd_int64 ((int64_t) e1, (int64_t) e2);
  e1 /= -d;
  e2 /= d;
  /* we will multiply row i1 by e2, and row i2 by e1 */

  t1 = 1;
  t2 = 1;
  t = 0;

  /* now perform the real merge */
  typerow_t *sum;
  sum = heap_alloc_row(i1, k1 + k2 - 1);

  int64_t e;
  while (t1 <= k1 && t2 <= k2) {
      if (r1[t1].id == r2[t2].id) {
	  /* as above, the exponent e below cannot overflow */
	  e = (int64_t) e2 * (int64_t) r1[t1].e + (int64_t) e1 * (int64_t) r2[t2].e;
	  if (e != 0) { /* exponents do not cancel */
	      ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
	      t++;
	      setCell(sum, t, r1[t1].id, e);
	    } else { // cancelation
	  	t1 ++;
	  	t2 ++;
	    }
	}
      else if (r1[t1].id < r2[t2].id)
	{
	  e = (int64_t) e2 * (int64_t) r1[t1].e;
	  ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
	  t++;
	  setCell(sum, t, r1[t1].id, e);
	  t1 ++;
	}
      else
	{
	  e = (int64_t) e1 * (int64_t) r2[t2].e;
	  ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
	  t++;
	  setCell(sum, t, r2[t2].id, e);
	  t2 ++;
	}
    }
  while (t1 <= k1) {
      e = (int64_t) e2 * (int64_t) r1[t1].e;
      ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
      t++;
      setCell(sum, t, r1[t1].id, e);
      t1 ++;
    }
  while (t2 <= k2) {
      e = (int64_t) e1 * (int64_t) r2[t2].e;
      ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
      t++;
      setCell(sum, t, r2[t2].id, e);
      // increase_weight (mat, r2[t2].id);
      t2 ++;
    }
  ASSERT(t <= k1 + k2 - 1);
  heap_resize_last_row(sum, t);
  heap_destroy_row(rows[i1]);
  rows[i1] = sum;
}
#endif

/* construct the D matrix (n' rows, n columns), where n is the number of rows
   of the original matrix M (output from purge), i.e., nrows. 
   Set eliminated_columns[j] == 1 if column j is eliminated */
static void
build_left_matrix(const char *outputname, const char *hisname, int bin)
{
	printf("Reading history file %s (2nd pass) and building left matrix\n", hisname);
	fflush(stdout);

	FILE * hisfile = fopen_maybe_compressed(hisname, "r");
	ASSERT_ALWAYS(hisfile != NULL);

	uint64_t addread = 0;
	char str[STRLENMAX];

	/* allocate identity matrix */
	typerow_t ** rowsL = malloc(nrows * sizeof(*rows));
	ASSERT_ALWAYS(rowsL != NULL);
	for (index_t i = 0; i < nrows; i++) {
		rowsL[i] = heap_alloc_row(i, 1);
		setCell(rowsL[i], 1, i, 1);
	}

	/* will print report at 2^10, 2^11, ... 2^23 computed primes and every
	 * 2^23 primes after that */
	stats_data_t stats;	/* struct for printing progress */
	stats_init(stats, stdout, &addread, 23, "Read", "row additions", "", "lines");

	nrows_small = nrows;
	int twomerge_mode = 1;
	uint64_t ntwomerge = 0;
	index_t left_nrows = nrows;
	while (fgets(str, STRLENMAX, hisfile)) {
		if (str[0] == '#')
			continue;

		addread++;

		if (stats_test_progress(stats))
			stats_print_progress(stats, addread, 0, 0, 0);

		if (str[strlen(str) - 1] != '\n') {
			fprintf(stderr, "Gasp: not a complete line!");
			fprintf(stderr, " I stop reading and go to the next phase\n");
			break;
		}

		index_t j;
		index_signed_t ind[MERGE_LEVEL_MAX], i0;
		int destroy;
		int ni = parse_hisfile_line(ind, str, &j);   // in sparse.c, mutualized with "normal" replay
	
		column_info[j] = 1;         // column has been eliminated

		if (ind[0] < 0) {
			destroy = 0;
			i0 = -ind[0] - 1;
		} else {
			destroy = 1;
			i0 = ind[0];
		}

		int twomerge = (ni == 2) && destroy;
		if (!twomerge)
			twomerge_mode = 0;

		if (twomerge_mode) {
			/* initial run of two merges: do them on R */
			add_row(rows, ind[1], i0, j);
			heap_destroy_row(rows[i0]);        // reclaim the memory
			rows[i0] = NULL;
			ntwomerge += 1;
			nrows_small -= 1;
		} else {
			/* normal case: history replayed on L */
			for (int k = 1; k < ni; k++)
				add_row(rowsL, ind[k], i0, j);
		}
		if (destroy) {
			heap_destroy_row(rowsL[i0]);        // reclaim the memory
			rowsL[i0] = NULL;
			left_nrows -= 1;
		}
	}
	stats_print_progress(stats, addread, 0, 0, 1);
	fclose_maybe_compressed(hisfile, hisname);
	printf("%" PRId64 " 2-merges done directly on R\n", ntwomerge);

	/* renumber columns of L to account for two-merges */
	index_t *renumber = malloc(nrows * sizeof(*renumber));
	ASSERT_ALWAYS(renumber != NULL);
	index_t acc = 0;
	for (index_t i = 0; i < nrows; i++) {
		if (rows[i] == NULL) {
			renumber[i] = UMAX(index_t);
		} else {
			renumber[i] = acc;
			acc += 1;
		}
	}
	ASSERT(acc == nrows_small);

	/* similar code below --- factorize? */
	for (index_t i = 0; i < nrows; i++) {            // renumber the columns
		if (rowsL[i] == NULL)
			continue;                        // row has been deleted
		for (index_t k = 1; k <= rowLength(rowsL, i); k++) {
			index_t j = rowCell(rowsL[i], k);
			ASSERT(column_info[j] != UMAX(index_t));
			#ifdef FOR_DL
				int32_t e = rowFullCell(rowsL[i], k).e; 
				setCell(rowsL[i], k, renumber[j], e);
			#else
				setCell(rowsL[i], k, renumber[j], 0);
			#endif
		}
	}
	free(renumber);

	/* output left matrix */
	flushSparse(outputname, rowsL, nrows, left_nrows, nrows_small, 0, bin);    // skip=0
	free(rowsL);
}

/******************************************************************************/
 
static void
build_right_matrix (const char *outputname, const char *idealsfilename, index_t skip, int bin)
{
	/* here: column_info[j] == 1   <====>   column has been eliminated
		 column_info[j] == 2   <====>   column is non-empty (not eliminated)
		 column_info[j] == 0   <====>   column is empty (not eliminated) */

	/* renumber the remaining non-empty columns (exclusive prefix-sum) */
	index_t sum = 0;
	for (uint64_t j = 0; j < ncols; j++) {
		if (column_info[j] != 2) {
			column_info[j] = UMAX(index_t);
			continue;
		}
		column_info[j] = sum;
		sum += 1;
	}
	printf("remaining (non-eliminated) columns : %" PRId64 "\n", (uint64_t) sum);

	/* 
	 * here: sum == number of non-eliminated columns.
	 *        column_info[j] == UMAX(...) ---> col j is out of the game
	 *        column_info[j] == k ---> col j becomes col k
	 */
	for (index_t i = 0; i < nrows; i++) {            // renumber the columns
		if (rows[i] == NULL)
			continue;                        // row has been deleted
		for (index_t k = 1; k <= rowLength(rows, i); k++) {
			index_t j = rowCell(rows[i], k);
			ASSERT(column_info[j] != UMAX(index_t));
			#ifdef FOR_DL
				int32_t e = rowFullCell(rows[i], k).e; 
				setCell(rows[i], k, column_info[j], e);
			#else
				setCell(rows[i], k, column_info[j], 0);
			#endif
		}
	}

	if (idealsfilename != NULL) {
		FILE *renumberfile = fopen_maybe_compressed (idealsfilename, "w");
		if (renumberfile == NULL) {
			fprintf (stderr, "Error while opening file to save permutation of ideals\n");
			exit(EXIT_FAILURE);
		}
		for (index_t j = 0; j < ncols; j++)
			if (column_info[j] == UMAX(index_t))
				fprintf(renumberfile, "# column %" PRIu64 " has been eliminated / is empty\n", (uint64_t) j);
			else
				fprintf (renumberfile, "%" PRIu64 " %" PRIx64 "\n",
					(uint64_t) column_info[j], (uint64_t) j);
		fclose(renumberfile);
	}

	/* output right matrix */
	flushSparse(outputname, rows, nrows, nrows_small, sum, skip, bin);
}

/******************************************************************************/

static void declare_usage(param_list pl)
{
	param_list_decl_usage(pl, "purged", "input purged file");
	param_list_decl_usage(pl, "his", "input history file");
	param_list_decl_usage(pl, "outL", "basename for left output matrices");
	param_list_decl_usage(pl, "outR", "basename for right output matrices");
#ifndef FOR_DL
	param_list_decl_usage(pl, "skip", "number of heaviest columns that go to the " "dense matrix (default " CADO_STRINGIZE(DEFAULT_MERGE_SKIP) ")");
#endif
	param_list_decl_usage(pl, "index", "file containing description of rows " "(relations-sets) of the matrix");
	param_list_decl_usage(pl, "ideals", "file containing correspondence between " "ideals and matrix columns");
	param_list_decl_usage(pl, "force-posix-threads", "force the use of posix threads, do not rely on platform memory semantics");
	param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
	param_list_decl_usage(pl, "ideals", "file containing correspondence between ideals and matrix columns");

	verbose_decl_usage(pl);
}

static void usage(param_list pl, char *argv0)
{
	param_list_print_usage(pl, argv0, stderr);
	exit(EXIT_FAILURE);
}

// We start from M_purged which is nrows x ncols;
int main(int argc, char *argv[])
{
	char *argv0 = argv[0];
        int skip = DEFAULT_MERGE_SKIP;
	double cpu0 = seconds();
	double wct0 = wct_seconds();
        int bin = -1;

#ifdef HAVE_MINGW
	_fmode = _O_BINARY;	/* Binary open for all files */
#endif

	setbuf(stdout, NULL);   // CB: where is setbuf? what's the purpose?
	setbuf(stderr, NULL);

#ifdef FOR_DL
	fprintf (stderr, "this is not ready. In particular, add_row must find the correct coefficients for linear combinations\n");
        exit (1);
#endif

	param_list pl;
	param_list_init(pl);
	declare_usage(pl);
	argv++, argc--;
	param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

	if (argc == 0)
		usage(pl, argv0);

	for (; argc;) {
		if (param_list_update_cmdline(pl, &argc, &argv)) {
			continue;
		}
		fprintf(stderr, "Unknown option: %s\n", argv[0]);
		usage(pl, argv0);
	}
	/* print command-line arguments */
	verbose_interpret_parameters(pl);
	param_list_print_command_line(stdout, pl);
	fflush(stdout);

	const char *purgedname = param_list_lookup_string(pl, "purged");
	const char *hisname = param_list_lookup_string(pl, "his");
	const char *sparseLname = param_list_lookup_string(pl, "outL");
	const char *sparseRname = param_list_lookup_string(pl, "outR");
	// const char *indexname = param_list_lookup_string(pl, "index");
	const char *idealsfilename = param_list_lookup_string(pl, "ideals");
	param_list_parse_int(pl, "skip", &skip);
	const char *path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");

	/* Some checks on command line arguments */
	if (param_list_warn_unused(pl)) {
		fprintf(stderr, "Error, unused parameters are given\n");
		usage(pl, argv0);
	}

	if (purgedname == NULL) {
		fprintf(stderr, "Error, missing -purged command line argument\n");
		usage(pl, argv0);
	}
	if (hisname == NULL) {
		fprintf(stderr, "Error, missing -his command line argument\n");
		usage(pl, argv0);
	}
	if (sparseLname == NULL || sparseRname == NULL) {
		fprintf(stderr, "Error, both -outL and --outR are required\n");
		usage(pl, argv0);
	}
#ifdef FOR_DL
	if (idealsfilename == NULL) {
		fprintf(stderr, "Error, missing -ideals command line argument\n");
		usage(pl, argv0);
	}
	ASSERT_ALWAYS(skip == 0);
#endif
        if (has_suffix(sparseLname, ".bin") || has_suffix(sparseLname, ".bin.gz"))
        {
          bin = 1;
          printf("# Output matrices will be written in binary format\n");
        }
        else
        {
          bin = 0;
          printf ("# Output matrices will be written in text format\n");
        }

	set_antebuffer_path(argv0, path_antebuffer);

	/* Read number of rows and cols on first line of purged file */
	uint64_t __nr, __nc;
	purgedfile_read_firstline(purgedname, &__nr, &__nc);  // uint64_t args
	if (__nr >= 4294967296UL) {
		fprintf(stderr, "Error, cannot handle 2^32 rows or more after purge\n");
		fprintf(stderr, "change ind_row from uint32_t to uint64_t in sparse.h\n");
		exit(EXIT_FAILURE);
	}
	printf("Purged matrix has %" PRIu64 " rows and %" PRIu64 " cols\n", __nr, __nc);
	fflush(stdout);
	nrows = __nr;
	ncols = __nc;

#if SIZEOF_INDEX == 4
	if (__nc >= UINT32_MAX) {
		fprintf(stderr, "You must recompile with -DSIZEOF_INDEX=8\n");
		exit(EXIT_FAILURE);
	}
#endif

	heap_setup();
	preread_history(hisname);

	/* load the relations in memory */
	read_purgedfile(purgedname);

	/* Read the history */
	printf("Building left matrix\n");
	build_left_matrix(sparseLname, hisname, bin);

	heap_reset();

	printf("Building right matrix\n");
	build_right_matrix(sparseRname, idealsfilename, skip, bin);

	printf("Cleaning up\n");
 	heap_clear();
	free(column_info);
	free(rows);
	param_list_clear(pl);
	print_timing_and_memory(stdout, cpu0, wct0);
	return 0;
}


// TODO : kill empty cols in purged mat
