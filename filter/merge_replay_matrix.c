#include "cado.h" // IWYU pragma: keep

#ifdef  __GNUC__
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef FOR_DL
#include "gcd.h"
#endif
#include "filter_config.h"
#include "merge_replay_matrix.h"
#include "merge_heap.h"
#include "sparse.h"
#include "macros.h"

/*****************************************************************************/

/* Initialize the sparse matrix mat. */
void
initMat (filter_matrix_t *mat, uint32_t skip)
{
  heap_setup(mat->heap);
  /* we start with cwmax = 2, and we increase it in mergeOneByOne() when
     the Markowitz queue is empty */
  mat->cwmax = 2;
  ASSERT_ALWAYS (mat->cwmax < 255); /* 255 is reserved for saturated values */
  mat->skip = skip;

  mat->p = NULL;
  mat->weight = 0;
  mat->tot_weight = 0;
  mat->rem_ncols = 0;
  mat->initial_ncols = mat->ncols;

  mat->rows = (typerow_t **) malloc (mat->nrows * sizeof (typerow_t *));
  ASSERT_ALWAYS (mat->rows != NULL);
  mat->wt = (col_weight_t *) malloc (mat->ncols * sizeof (col_weight_t));
  ASSERT_ALWAYS (mat->wt != NULL);
  memset (mat->wt, 0, mat->ncols * sizeof(col_weight_t));
  mat->p = NULL; /* recompress() assumes mat->p = 0 at the beginning,
		    in which case it just renumbers */
}

void
clearMat (filter_matrix_t *mat)
{
  free (mat->rows);
  free (mat->wt);
  heap_clear(mat->heap);
}

int cmp_u64(const uint64_t * a, const uint64_t * b)
{
    return (*a > *b) - (*b > *a);
}

void
print_row(filter_matrix_t *mat, index_t i)
{
    fprintRow (stdout, mat->rows[i]);
}
