#ifndef MERGE_HEAP_H_
#define MERGE_HEAP_H_

#include "omp_proxy.h"	        // omp_test_lock and the likes
#include "merge_types.h"        // typerow_t, index_t

/* 
 * Custom memory management for arrays of typerow_t (i.e. sparse matrices)
 * where rows are frequently allocated/deleted.
 */


#ifdef __cplusplus
extern "C" {
#endif

#define PAGE_DATA_SIZE ((1<<18) - 4) /* seems to be optimal for RSA-512 */

struct page_t {
        struct pagelist_t *list;     /* the pagelist_t structure associated with this page */
        int i;                       /* page number, for debugging purposes */
        int generation;              /* pass in which this page was filled. */
        int ptr;                     /* data[ptr:PAGE_DATA_SIZE] is available*/
        typerow_t data[PAGE_DATA_SIZE];
};

// linked list of pages (doubly-linked for the full pages, simply-linked for the empty pages)
struct pagelist_t {
        struct pagelist_t *next;
        struct pagelist_t *prev;
        struct page_t *page;
};

struct __heap_ctx_struct {
	int canary;
	int n_pages, n_full_pages, n_empty_pages;
	struct pagelist_t headnode;                  // dummy node for the list of full pages
	struct pagelist_t *full_pages, *empty_pages; // head of the page linked lists
	struct page_t **active_page; /* active page of each thread */
	long long *heap_waste;       /* space wasted, per-thread. Can be negative! The sum over all threads is correct. */
	int current_generation;
	omp_lock_t lock;
};

typedef struct __heap_ctx_struct heapctx_t[1];

/* Setup data structures. Must be called outside of a parallel region */
void heap_setup(heapctx_t ctx);

int heap_config_get_PAGE_DATA_SIZE();

/* Allocate space for row i, holding s coefficients in row[1:s+1] (row[0] == s). 
	This writes s in row[0].  The size of the row must not change afterwards.
	Thread-safe.  Usually very fast. */
typerow_t * heap_alloc_row (heapctx_t ctx, index_t i, size_t s);

/* Shrinks the last allocated row.  It must be the last one allocated by this
	thread. Lock-free and constant-time. */
void heap_resize_last_row (heapctx_t ctx, typerow_t *row, index_t new_size);

/* Given the pointer provided by heap_alloc_row, mark the row as deleted.
	The memory is not released.  Lock-free and constant-time. */
void heap_destroy_row (heapctx_t ctx, typerow_t *row);

/* Reclaim all the memory occupied by destroyed rows. 
	Must be called outside of a parallel region.  May move rows in memory and update row pointers. */
void heap_garbage_collection(heapctx_t ctx, typerow_t **rows);

/* Deallocate all rows. Does not release memory (it remains available for row allocation) */
void heap_reset (heapctx_t ctx);

/* release all memory. This is technically not necessary, because the "malloc"
	allocations are internal to the process, and all space allocated to the
	process is reclaimed by the OS on termination. However, doing this enables
	valgrind to check the absence of leaks.
*/
void heap_clear (heapctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif	/* MERGE_HEAP_H_ */
