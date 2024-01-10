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
// #include <libgen.h>
#include "purgedfile.h"		// for purgedfile_read_firstline
#include "typedefs.h"		// for index_t, ideal_merge_t, index_signed_t
#include "filter_config.h"
#include "filter_io.h"		// earlyparsed_relation_ptr
#include "fix-endianness.h"	// fwrite32_little
#include "gzip.h"		// fopen_maybe_compressed
#include "misc.h"		// derived_filename, has_suffix
#include "params.h"		// param_list_parse_*
#include "stats.h"		// stats_data_t
#include "timing.h"		// seconds
#include "verbose.h"		// verbose_decl_usage
#include "portability.h"	// strdup  // IWYU pragma: keep
#include "macros.h"
#include "merge_heap.h"
#include "gcd.h"
#include "merge_replay_matrix.h" // rowCell

#define DEBUG 0

/* purged matrix */
heapctx_t heap;
typerow_t **rows;       // purged matrix
index_t nrows;          // #rows for the purged matrix
index_t ncols;          // #cols for the purged matrix
index_t nrows_small;    // #rows for the purged matrix once 2-merges are performed        
index_t left_nrows;     // #rows for the left (& product) matrix

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
        if (!bin) 
        	fprintf(smatfile, "%" PRIu64 " %" PRIu64 "\n",
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
			for (index_t k = 1; k <= rowLength(rows, i); k++) {
				index_t j = rowCell(rows[i], k);
				assert(j < ncols);
				if (j < skip) {
					dw++;
					DW++;
				} else {
					sw++;
					W++;
				}
			}
			if (bin) {
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

/*************************************************************************/

/* COPIED-PASTED from merge.c. REFACTORING PLAN: take rows[] and weights[] as
   arugments, update weights only if not NULL */


#ifndef FOR_DL
/* special code for factorization */

// row[i1] <--- row[i1] + row[i2]
static void
add_row(typerow_t **rows, index_t i1, index_t i2)
{
	typerow_t *r1 = rows[i1];
  	typerow_t *r2 = rows[i2];
	index_t k1 = rowLength(rows, i1);
	index_t k2 = rowLength(rows, i2);
	index_t t1 = 1;            // index in rows[i1]
	index_t t2 = 1;            // index in rows[i2]
	index_t t = 0;             // index in the sum

	/* fast-track : don't precompute the size */
	typerow_t *sum = heap_alloc_row(heap, i1, k1 + k2);

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
	heap_resize_last_row(heap, sum, t);
	heap_destroy_row(heap, r1);
	rows[i1] = sum;
	return;
}
#else /* FOR_DL */

#define INT32_MIN_64 (int64_t) INT32_MIN
#define INT32_MAX_64 (int64_t) INT32_MAX

// row[i1] <--- row[i1] * e2 + row[i2] * e1
static void
add_row(typerow_t **rows, index_t i1, int32_t e2, index_t i2, int32_t e1)
{
 	typerow_t *r1 = rows[i1];
  	typerow_t *r2 = rows[i2];
 	index_t k1 = rowLength(rows, i1);
	index_t k2 = rowLength(rows, i2);
	index_t t1 = 1;            // index in rows[i1]
	index_t t2 = 1;            // index in rows[i2]
	index_t t = 0;             // index in the sum

  /* we always check that e1 and e2 are not zero, in order to prevent from zero
     exponents that would come from exponent overflows in previous merges */
	ASSERT_ALWAYS (e1 != 0 && e2 != 0);

  t1 = 1;
  t2 = 1;
  t = 0;

  /* now perform the real merge */
  typerow_t *sum;
  sum = heap_alloc_row(heap, i1, k1 + k2);

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
      t2 ++;
    }
  ASSERT(t <= k1 + k2);
  heap_resize_last_row(heap, sum, t);
  heap_destroy_row(heap, rows[i1]);
  rows[i1] = sum;
}
#endif

#define STRLENMAX 2048

/******************************************************************************/

/* similar (but not identical) to the same function in replay.c */
static void
writeIndex(const char *indexname, typerow_t **rows, index_t small_nrows)
{
    FILE *indexfile = NULL;
    indexfile = fopen_maybe_compressed(indexname, "w");
    ASSERT_ALWAYS (indexfile != NULL);
    printf("Writing index to %s\n", indexname);

    fprintf(indexfile, "%" PRIu64 "\n", (uint64_t) small_nrows);

    for (index_t i = 0; i < small_nrows; ++i) {
    	typerow_t MAYBE_UNUSED *row = rows[i];
    	index_t row_length = rowLength(rows, i);
        ASSERT (row_length > 0);
        fprintf(indexfile, "%" PRIu64, (uint64_t) row_length);
        for (unsigned int j = 1; j <= row_length; j++) {
#ifdef FOR_DL
            fprintf(indexfile, " %" PRIx64 ":%d",
                    (uint64_t) row[j].id, row[j].e);
#else
            fprintf(indexfile, " %" PRIx64 "", (uint64_t) row[j]);
#endif
        }
        fprintf(indexfile, "\n");
    }
    fclose_maybe_compressed(indexfile, indexname);
}


/*  Set column_info[j] == { 1 if column j is eliminated 
 *                        { 0  otherwise
 *  Output the "index" file --- a more verbose version of the left matrix
 */
static void
preread_history(const char *hisname, const char *indexname)
{
	printf("Reading history file %s (1st pass) and producing index\n", hisname);
	fflush(stdout);

	FILE * hisfile = fopen_maybe_compressed(hisname, "r");
	ASSERT_ALWAYS(hisfile != NULL);

	/* allocate auxiliary data (per column) */
	index_t n_elim = 0;
	column_info = malloc(sizeof(*column_info) * ncols);
	ASSERT_ALWAYS(column_info != NULL);
	for (index_t i = 0; i < ncols; i++)
		column_info[i] = 0;

	/* allocate identity matrix for index */
	typerow_t **rows_index = malloc(nrows * sizeof(*rows));
	ASSERT_ALWAYS(rows_index != NULL);
	for (index_t i = 0; i < nrows; i++) {
		rows_index[i] = heap_alloc_row(heap, i, 1);
		setCell(rows_index[i], 1, i, 1);
	}
	/** BEGIN NOT DRY (w.r.t. build left matrix) ***/

	/* will print report at 2^10, 2^11, ... 2^23 computed primes and every
	 * 2^23 primes after that */
	index_t left_nrows = nrows;
	uint64_t addread = 0;
	char str[STRLENMAX];
	stats_data_t stats;	/* struct for printing progress */
	stats_init(stats, stdout, &addread, 23, "Read", "row additions", "", "lines");
	while (fgets(str, STRLENMAX, hisfile)) {
		addread++;
		if (stats_test_progress(stats))
			stats_print_progress(stats, addread, 0, 0, 0);
		if (str[strlen(str) - 1] != '\n') {
			fprintf(stderr, "Gasp: not a complete line!");
			fprintf(stderr, " I stop reading and go to the next phase\n");
			break;
		}

        	uint64_t i1, i2, j;
#ifdef FOR_DL
        	int32_t e1, e2;
#endif
        	int k;
        	switch(str[0]) {
        	case '#':
        	case '!':
        		break;
        	case '|':
			k = sscanf(str + 1, "%" SCNd64, &j);
        		ASSERT(k == 1);
        		/* mark column j as being eliminated */
			ASSERT(column_info[j] != 1);
			column_info[j] = 1;
			n_elim += 1;
        	        break;
        	case '-':        // destroy row
        	        k = sscanf(str + 1, "%" SCNd64, &i1);
        	        ASSERT(k == 1);
        	        heap_destroy_row(heap, rows_index[i1]);
			left_nrows -= 1;
        	        rows_index[i1] = NULL;
        	        break;
#ifndef FOR_DL
        	case '+':        // row addition
                	k = sscanf(str + 1, "%" SCNd64 " %" SCNd64, &i1, &i2);
                	ASSERT(k == 2);
                	add_row(rows_index, i1, i2);
                	break;
        	case '*':
        	        fprintf (stderr, "coefficients not allowed in factorization mode\n");
            		exit (EXIT_FAILURE);	
#else   // FOR_DL
            	case '+':
            		fprintf (stderr, "absence of coefficients not allowed in DLP mode\n");
            		exit (EXIT_FAILURE);	
            	case '*':       // row addition with multiplicative coefficients
        	        k = sscanf(str + 1, "%" SCNd64 " %" PRId32 " %" SCNd64 " %" PRId32, &i1, &e1, &i2, &e2);
                	ASSERT(k == 4);
                	add_row(rows_index, i1, e1, i2, e2);
        	        break;
#endif
        	default:
			fprintf (stderr, "parse error in history file\n");
            		exit (EXIT_FAILURE);
        	}
	}
	stats_print_progress(stats, addread, 0, 0, 1);
	fclose_maybe_compressed(hisfile, hisname);

	/** END NOT DRY (w.r.t. build left matrix) ***/
	printf("%" PRId64 " eliminated columns\n", (uint64_t) n_elim);

	/* emit index */
	index_t j = 0;
	for (index_t i = 0; i < nrows; i++)  /* stack non-empty rows in index */
		if (rows_index[i] != NULL) {
			rows_index[j] = rows_index[i];
			j += 1;
		}
	ASSERT_ALWAYS(j == left_nrows);
	writeIndex(indexname, rows_index, left_nrows);
	free(rows_index);
}


/*
 * The following code is very similar to replay.c
 * It loads the whole matrix in memory.
 */
void * read_purged_row (void MAYBE_UNUSED *context_data, earlyparsed_relation_ptr rel)
{
	/* 
	 * 1st pass, set scratch (parity of #occurence of each column)
	 * This is only useful in factoring
	 */
	#ifndef FOR_DL
		for (unsigned int j = 0; j < rel->nb; j++) {
			index_t h = rel->primes[j].h;
			if (column_info[h] == 1)          // column was eliminated
				continue;
			scratch[h] ^= 1;
		}
	#endif

	/* 2nd pass, copy row and reset scratch */
	typerow_t buf[UMAX(weight_t)];
	unsigned int nb = 0;
	for (unsigned int j = 0; j < rel->nb; j++) {
		index_t h = rel->primes[j].h;
		if (column_info[h] == 1)          // column was eliminated
			continue;
		#ifndef FOR_DL 
			/* factoring: test parity of #occurence of this ideal */
			if (scratch[h] == 0)
				continue; // ideal appeared an even number of times
		#endif
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

	rows[rel->num] = heap_alloc_row(heap, rel->num, nb);  // in merge_heap.c
	compressRow (rows[rel->num], buf, nb);                // in sparse.c
  	return NULL;
}

/* update rows, column_info. 
 * set column_info[j] == 1   <====>   column has been eliminated
 *     column_info[j] == 2   <====>   column is non-empty (not eliminated)
 *     column_info[j] == 0   <====>   column is empty (not eliminated) 
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

	uint64_t empty_cols = 0;
	for (index_t j = 0; j < ncols; j++)
		if (column_info[j] == 0)
			empty_cols += 1;
	printf("Found %" PRId64 " empty columns in the purged matrix\n", empty_cols);
}

/* construct the left matrix (n' rows, n columns), where n is the number of rows
   of the original matrix M (output from purge), i.e., nrows. 
   This is similar to the doAllAdds() function in replay.c
*/
static void
build_left_matrix(const char *outputLname, const char *outputPname, const char *hisname, index_t skip, int bin)
{
	printf("Reading history file %s (2nd pass), building left matrix\n", hisname);
	fflush(stdout);

	FILE * hisfile = fopen_maybe_compressed(hisname, "r");
	ASSERT_ALWAYS(hisfile != NULL);

	/* allocate identity matrix for L */
	typerow_t ** rowsL = malloc(nrows * sizeof(*rowsL));
	ASSERT_ALWAYS(rowsL != NULL);
	for (index_t i = 0; i < nrows; i++) {
		rowsL[i] = heap_alloc_row(heap, i, 1);
		setCell(rowsL[i], 1, i, 1);
	}

	/* allocate and fill "product matrix" for the skipped (=dense) columns */
	typerow_t ** rowsP = malloc(nrows * sizeof(*rowsP));
	ASSERT_ALWAYS(rowsP != NULL);
	for (index_t i = 0; i < nrows; i++) {
		int l = 0;  /* #skipped entries */
		for (index_t k = 1; k <= rowLength(rows, i); k++)
			if (rowCell(rows[i], k) < skip)
				l += 1;
		rowsP[i] = heap_alloc_row(heap, i, l);
		l = 1;
		for (index_t k = 1; k <= rowLength(rows, i); k++)
			if (rowCell(rows[i], k) < skip) {
				rowsP[i][l] = rows[i][k];
				l += 1;
			}
	}

	/* read and process history */
	/* will print report at 2^10, 2^11, ... 2^23 computed primes and every
	 * 2^23 primes after that */
	uint64_t ntwomerge = 0;
	nrows_small = nrows;
	int twomerge_mode = 1;
	left_nrows = nrows;
	uint64_t addread = 0;
	char str[STRLENMAX];	stats_data_t stats;	/* struct for printing progress */
	stats_init(stats, stdout, &addread, 23, "Read", "row additions", "", "lines");
	
	while (fgets(str, STRLENMAX, hisfile)) {
		addread++;
		if (stats_test_progress(stats))
			stats_print_progress(stats, addread, 0, 0, 0);
		/* sanity check */
		if (str[strlen(str) - 1] != '\n') {
			fprintf(stderr, "Gasp: not a complete line!");
			fprintf(stderr, " I stop reading and go to the next phase\n");
			break;
		}

		// TODO: garbage-collect. Distinguish the different heaps

        	uint64_t i1, i2;
#ifdef FOR_DL
        	int32_t e1, e2;
#endif
        	int k;
        	switch(str[0]) {
        	case '#':
        		break;
        	case '|':
        		ntwomerge += 1;
        	        break;
        	case '!':
        		twomerge_mode = 0;
        		break;
        	case '-':        // destroy row
			k = sscanf(str + 1, "%" SCNd64, &i1);
			ASSERT(k == 1);
			if (twomerge_mode) {
	       	        	heap_destroy_row(heap, rows[i1]);
             			rows[i1] = NULL;
				nrows_small -= 1;
			}
			heap_destroy_row(heap, rowsL[i1]);
			heap_destroy_row(heap, rowsP[i1]);
			rowsL[i1] = NULL;
			left_nrows -= 1;
			break;
#ifndef FOR_DL
        	case '+':        // row addition
                	k = sscanf(str + 1, "%" SCNd64 " %" SCNd64, &i1, &i2);
                	ASSERT(k == 2);
                	if (twomerge_mode)
                		add_row(rows, i1, i2);
                	else
                		add_row(rowsL, i1, i2);
                	add_row(rowsP, i1, i2);
                	break;
        	case '*':
        	        fprintf (stderr, "coefficients not allowed in factorization mode\n");
            		exit (EXIT_FAILURE);	
#else   // FOR_DL
            	case '+':
            		fprintf (stderr, "absence of coefficients not allowed in DLP mode\n");
            		exit (EXIT_FAILURE);	
            	case '*':       // row addition with multiplicative coefficients
        	        k = sscanf(str + 1, "%" SCNd64 " %" PRId32 " %" SCNd64 " %" PRId32, &i1, &e1, &i2, &e2);
                	ASSERT(k == 4);
                	if (twomerge_mode)
                		add_row(rows, i1, e1, i2, e2);
			else
                		add_row(rowsL, i1, e1, i2, e2);
                	add_row(rowsP, i1, e1, i2, e2);
        	        break;
#endif
        	default:
			fprintf (stderr, "parse error in history file\n");
            		exit (EXIT_FAILURE);
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
			ASSERT(renumber[j] != UMAX(index_t));
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
	flushSparse(outputLname, rowsL, nrows, left_nrows, nrows_small, 0, bin);    // skip=0
	free(rowsL);

	/* output product matrix */
	flushSparse(outputPname, rowsP, nrows, left_nrows, skip, skip, bin);
	free(rowsP);
}

/******************************************************************************/
 
static void
export_right_matrix (const char *outputname, const char *idealsname, index_t skip, int bin)
{
	/* here: column_info[j] == 1   <====>   column has been eliminated
	 *       column_info[j] == 0   <====>   column is empty (not eliminated) 
	 *       column_info[j] == 2   <====>   column is non-empty (not eliminated) 
	 */

	/* renumber the remaining non-empty columns (exclusive scan) */
	index_t sum = 0;
	for (uint64_t j = 0; j < ncols; j++) {
		if (column_info[j] != 2) {
			column_info[j] = UMAX(index_t);
			continue;
		}
		column_info[j] = sum;
		sum += 1;
	}
	printf("remaining (non-eliminated, non-empty) columns : %" PRId64 "\n", (uint64_t) sum);

	/* 
	 * here: sum == number of non-eliminated columns.
	 *        column_info[j] == UMAX(...) ---> col j is out of the game
	 *        column_info[j] == k         ---> col j becomes col k
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

	if (idealsname != NULL) {
		FILE *renumberfile = fopen_maybe_compressed (idealsname, "w");
		if (renumberfile == NULL) {
			fprintf (stderr, "Error while opening file to save permutation of ideals\n");
			exit(EXIT_FAILURE);
		}
		printf("Saving ideal renumbering in %s\n", idealsname);
		fprintf (renumberfile, "# %" PRIu64 "\n", (uint64_t) sum);
		for (index_t j = 0; j < ncols; j++)
			if (column_info[j] != UMAX(index_t))
			// 	fprintf(renumberfile, "# column %" PRIu64 " has been eliminated / is empty\n", (uint64_t) j);
			// else
				fprintf (renumberfile, "%" PRIu64 " %" PRIx64 "\n",
					(uint64_t) column_info[j], (uint64_t) j);
		fclose_maybe_compressed(renumberfile, idealsname);
	}

	/* output right matrix */
	flushSparse(outputname, rows, nrows, nrows_small, sum, skip, bin);
}

/******************************************************************************/

static void declare_usage(param_list pl)
{
	param_list_decl_usage(pl, "purged", "input purged file");
	param_list_decl_usage(pl, "his", "input history file");
	param_list_decl_usage(pl, "out", "basename for output matrices");
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

/* convert dirname/basename into dirname/Xbasename */
static void
make_out (char *s, const char *out, const char X)
{
	char *tmp;
	switch (X) {
	case 'L':
  		tmp = derived_filename(out, "L", ".sparse.bin");
		break;
	case 'R':
  		tmp = derived_filename(out, "R", ".sparse.bin");
  		break;
  	default:
  		ASSERT(0);
  	}
  	strcpy(s, tmp);
  	free(tmp);
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

#if 0
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
	const char *indexname = param_list_lookup_string(pl, "index");
	const char *sparsename = param_list_lookup_string(pl, "out");
        char sparseLname[1024], sparseRname[1024];
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
	if (sparsename == NULL) {
		fprintf(stderr, "Error, -out is required\n");
		usage(pl, argv0);
	}
	if (indexname == NULL) {
		fprintf(stderr, "Error, --index is required\n");
		usage(pl, argv0);
	}
#ifdef FOR_DL
	if (idealsfilename == NULL) {
		fprintf(stderr, "Error, missing -ideals command line argument\n");
		usage(pl, argv0);
	}
	ASSERT_ALWAYS(skip == 0);
#endif
        if (has_suffix(sparsename, ".bin") || has_suffix(sparsename, ".bin.gz")) {
          bin = 1;
          printf("# Output matrices will be written in binary format\n");
        } else {
          bin = 0;
          printf ("# Output matrices will be written in text format\n");
        }
	make_out(sparseLname, sparsename, 'L');
	make_out(sparseRname, sparsename, 'R');

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

	/* read history ; mark eliminated columns ; output index */
	heap_setup(heap);
	preread_history(hisname, indexname);

	/* load the relations in memory */
	heap_reset(heap);
	read_purgedfile(purgedname);

	printf("Building left & (dense) product matrix\n");
	build_left_matrix(sparseLname, sparsename, hisname, skip, bin);

	printf("Building right matrix\n");
	export_right_matrix(sparseRname, idealsfilename, skip, bin);

	printf("Cleaning up\n");
 	heap_clear(heap);
	free(column_info);
	free(rows);
	param_list_clear(pl);
	print_timing_and_memory(stdout, cpu0, wct0);
	return 0;
}


// TODO : kill empty cols in purged mat
