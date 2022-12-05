/* replay --- replaying history of merges to build the sparse matriCES in
              Kleinjung's "double matrix" idea.

Copyright 2008-2019 Francois Morain, Emmanuel Thome, Paul Zimmermann,
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
#include "gcd.h"                // gcd_int64

#define DEBUG 0

/* Save a sparse matrix to the filesystem. Write column and row weights in
   companion files. If skip > 0, also create a "dense" part in the same format. */ 
static unsigned long flushSparse(const char *sparsename, typerow_t ** rows, 
        index_t nrows, index_t ncols, index_t skip)
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
	} suffixes = {
		 .ext = ".bin",
		 .smat = "bin",
		 .srw = "rw.bin",
		 .scw = "cw.bin",
		 .dmat = "dense.bin",
		 .drw = "dense.rw.bin",
		 .dcw = "dense.cw.bin",
	};
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
	if (has_suffix(base, suffixes.ext)) { /* strip suffix if given */
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

	smatname = derived_filename(base, suffixes.smat, zip);
	srwname = derived_filename(base, suffixes.srw, zip);
	scwname = derived_filename(base, suffixes.scw, zip);
	smatfile = fopen(smatname, wmode);
	srwfile = fopen(srwname, wmode);

	if (skip) {
		/* arrange so that we don't get file names like .sparse.dense */
		char *dbase = strdup(base);
		char *tmp = strstr(dbase, ".sparse");
		if (tmp)
			memmove(tmp, tmp + 7, strlen(tmp + 7) + 1);
		dmatname = derived_filename(dbase, suffixes.dmat, zip);
		drwname = derived_filename(dbase, suffixes.drw, zip);
		dcwname = derived_filename(dbase, suffixes.dcw, zip);
		dmatfile = fopen(dmatname, wmode);
		drwfile = fopen(drwname, wmode);
		free(dbase);
	}

	for (index_t i = 0; i < nrows; i++) {
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
			fwrite32_little(&sw, 1, smatfile);
			if (srwfile)
				fwrite32_little(&sw, 1, srwfile);
			if (skip)
				fwrite32_little(&dw, 1, dmatfile);
			if (skip)
				fwrite32_little(&dw, 1, drwfile);

			for (index_t j = 1; j <= rowLength(rows, i); j++) {
				ASSERT_ALWAYS(rowCell(rows[i], j) <= (index_t) UINT32_MAX);
				uint32_t x = rowCell(rows[i], j);
				if (srwfile)
					weights[x]++;
				if (x < skip) {
					ASSERT_ALWAYS(skip);
				fwrite32_little(&x, 1, dmatfile);
				} else {
					x -= skip;
					fwrite32_little(&x, 1, smatfile);
#ifdef FOR_DL
					/* exponents are always int32_t */
					uint32_t e = (uint32_t) rows[i][j].e;
					fwrite32_little(&e, 1, smatfile);
#endif
					
				}
			}
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

typedef struct {
	typerow_t **mat;
	index_t ncols;
	index_t col0;
	index_t colmax;
} replay_read_data_t;


void *fill_in_rows(void *context_data, earlyparsed_relation_ptr rel)
{
	replay_read_data_t *data = (replay_read_data_t *) context_data;
	typerow_t buf[UMAX(weight_t)];

	unsigned int nb = 0;
	for (unsigned int j = 0; j < rel->nb; j++) {
		index_t h = rel->primes[j].h;
		if (h < data->col0 || h >= data->colmax)
			continue;
		nb++;
#ifdef FOR_DL
		exponent_t e = rel->primes[j].e;
		buf[nb] = (ideal_merge_t) {
		.id = h,.e = e};
#else
		ASSERT_ALWAYS(rel->primes[j].e == 1);
		buf[nb] = h;
#endif
		ASSERT(h < data->ncols);
	}
#ifdef FOR_DL
	buf[0].id = nb;
#else
	buf[0] = nb;
#endif

	qsort(&(buf[1]), nb, sizeof(typerow_t), cmp_typerow_t);

	data->mat[rel->num] = mallocRow(nb + 1);
	compressRow(data->mat[rel->num], buf, nb);

	return NULL;
}

/* if for_msieve=1, generate the *.cyc file needed by msieve to construct
   its matrix, which is of the following (binary) format:
      small_nrows
      n1 i1 i2 ... in1
      n2 j1 j2 ... jn2
      ...
      nk ...
   where each value is stored as a 32-bit integer (no linebreak),
   small_nrows is the number of relation-sets of the matrix,
   n1 is the number of relations in the first relation-set,
   i1 is the index of the first relation in the first relation-set
   (should correspond to line i1+2 in *.purged.gz), and so on */

MAYBE_UNUSED
static void read_purgedfile(typerow_t ** mat, const char *filename, index_t nrows, index_t ncols, index_t col0, index_t colmax)
{
	index_t nread;
	printf("Reading sparse matrix from %s\n", filename);
	fflush(stdout);
	char *fic[2] = { (char *)filename, NULL };
	replay_read_data_t tmp = (replay_read_data_t) {
		.mat = mat,
		.ncols = ncols,
		.col0 = col0,
		.colmax = colmax,
	};
	nread = filter_rels(fic, (filter_rels_callback_t) & fill_in_rows, &tmp, EARLYPARSE_NEED_INDEX, NULL, NULL);
	ASSERT_ALWAYS(nread == nrows);
}


// i1 += i2
// j is the index of the column that is used for
// pivoting in the case of DL. Then, the operation is
//   i1 = e2*i1 + e1*i2
// where e1 and e2 are adjusted so that the j-th column is zero in i1.

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
		sum[++t] = r2[t2];
		t1 += 1;
	}
	ASSERT(t <= k1 + k2 - 1);
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


static void
build_left_matrix(const char *matrixname, const char *hisname, index_t nrows, 
	index_t ncols, int skip, index_t Nmax)
{
	FILE * hisfile = fopen_maybe_compressed(hisname, "r");
	ASSERT_ALWAYS(hisfile != NULL);

	/* read merges in the *.merge.his file and replay them */
	uint64_t addread = 0;
	char str[STRLENMAX];

	printf("Reading row additions\n");
	fflush(stdout);

	/* allocate identity matrix */
	typerow_t ** rows = malloc(nrows * sizeof(*rows));
	ASSERT_ALWAYS(rows != NULL);
	for (index_t i = 0; i < nrows; i++) {
		rows[i] = heap_alloc_row(i, 1);
		rows[i][1] = i;
	}

	/* will print report at 2^10, 2^11, ... 2^23 computed primes and every
	 * 2^23 primes after that */
	stats_data_t stats;	/* struct for printing progress */
	stats_init(stats, stdout, &addread, 23, "Read", "row additions", "", "lines");
	
	int left_nrows = nrows;
	while (fgets(str, STRLENMAX, hisfile) && nrows >= Nmax) {
		if (str[0] == '#')
			continue;

		addread++;

		if (stats_test_progress(stats))
			stats_print_progress(stats, addread, 0, 0, 0);

		if (str[strlen(str) - 1] != '\n') {
			fprintf(stderr, "Gasp: not a complete a line!");
			fprintf(stderr, " I stop reading and go to the next phase\n");
			break;
		}
		
		index_t j;
		index_signed_t ind[MERGE_LEVEL_MAX], i0;
		int destroy;
		int ni = parse_hisfile_line(ind, str, &j);   // in sparse.c, mutualized with "normal" replay
		if (ind[0] < 0) {
			destroy = 0;
			i0 = -ind[0] - 1;
		} else {
			destroy = 1;
			i0 = ind[0];
		}

		for (int k = 1; k < ni; k++)
			add_row(rows, ind[k], i0, j);

		if (destroy) {
			heap_destroy_row(rows[i0]);
			rows[i0] = NULL;
			left_nrows -= 1;
		}			
	}
	stats_print_progress(stats, addread, 0, 0, 1);
	fclose_maybe_compressed(hisfile, hisname);

	/* output left matrix */
	flushSparse(matrixname, rows, nrows, ncols, skip);
}



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
	param_list_decl_usage(pl, "Nmax", "stop at Nmax number of rows (default 0)");
#ifndef FOR_DL
	param_list_decl_usage(pl, "col0", "print only columns with index >= col0");
	param_list_decl_usage(pl, "colmax", "print only columns with index < colmax");
#endif
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
	uint64_t Nmax = 0;
	uint64_t nrows, ncols;
        int skip = DEFAULT_MERGE_SKIP;
	double cpu0 = seconds();
	double wct0 = wct_seconds();

#ifdef HAVE_MINGW
	_fmode = _O_BINARY;	/* Binary open for all files */
#endif

	setbuf(stdout, NULL);   // CB: where is setbuf? what's the purpose?
	setbuf(stderr, NULL);

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
	// const char *idealsfilename = param_list_lookup_string(pl, "ideals");
#ifndef FOR_DL
	param_list_parse_int(pl, "skip", &skip);
#endif
	param_list_parse_uint64(pl, "Nmax", &Nmax);
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
	if (sparseLname == NULL && sparseRname == NULL) {
		fprintf(stderr, "Error, at least one of -outL or --outR is required\n");
		usage(pl, argv0);
	}
#ifdef FOR_DL
	if (idealsfilename == NULL) {
		fprintf(stderr, "Error, missing -ideals command line argument\n");
		usage(pl, argv0);
	}
	ASSERT_ALWAYS(skip == 0);
#endif
	printf("# Output matrices will be written in binary format\n");

	set_antebuffer_path(argv0, path_antebuffer);

	/* Read number of rows and cols on first line of purged file */
	purgedfile_read_firstline(purgedname, &nrows, &ncols);
	if (nrows >= 4294967296UL) {
		fprintf(stderr, "Error, cannot handle 2^32 rows or more after purge\n");
		fprintf(stderr, "change ind_row from uint32_t to uint64_t in sparse.h\n");
		exit(EXIT_FAILURE);
	}
	printf("Sparse matrix has %" PRIu64 " rows and %" PRIu64 " cols\n", nrows, ncols);
	fflush(stdout);

#if SIZEOF_INDEX == 4
	if (ncols >= UINT32_MAX) {
		fprintf(stderr, "You must recompile with -DSIZEOF_INDEX=8\n");
		exit(EXIT_FAILURE);
	}
#endif

	heap_setup();

	/* Read the matrix from purgedfile */
	if (sparseLname != NULL) {
		printf("Building left matrix\n");
		build_left_matrix(sparseLname, hisname, nrows, ncols, skip, Nmax);
	}

	/* Read the matrix from purgedfile */
	if (sparseRname != NULL) {
		// WIP
		// read_purgedfile(rows, purgedname, nrows, ncols, col0, colmax);
		printf("The biggest index appearing in a relation is %" PRIu64 "\n", ncols);
		fflush(stdout);
	}

	// fasterVersion(newrows, sparseRname, indexname, hisname, nrows, ncols, skip, idealsfilename, Nmax);

	param_list_clear(pl);
	print_timing_and_memory(stdout, cpu0, wct0);
	return 0;
}
