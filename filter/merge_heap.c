#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"
#include "merge_heap.h"
#include "merge_replay_matrix.h"  // typerow_t, index_t, 
#include "sparse.h" // rowCell
#include "omp_proxy.h"
#include "timing.h"  // seconds

/*************************** heap structures *********************************/

/* 
  Threads allocate PAGES of memory of a fixed size to store rows, and each
  thread has a single ACTIVE page in which it writes new rows.  When the
  active page is FULL, the thread grabs an EMPTY page that becomes its new
  active page.  

  In each page, there is a pointer [[ptr]] to the beginning of the free space.
  To allocate [[b]] bytes for a new row, it suffices to note the current
  value of[[ptr]] and then to increase it by [[b]] --- if this would overflow
  the current page, then it is marked as FULL and a new active page is
  obtained.  

  Rows are stored along with their number, their size and the list of their
  coefficients.  To delete a row, we just mark it as deleted by setting its
  number to -1.  Thus, row allocation and deallocation are thread-local
  operations that are very fast.  The last allocated row can easily be
  shrunk (by diminishing [[ptr]]).

  After each merge pass, memory is garbage-collected.  All threads do the
  following procedure in parallel, while possible: grab a FULL page that was
  not created during this pass (because there is no free space in it); copy 
  all non-deleted rows to the current ACTIVE page; mark the old FULL page as 
  EMPTY.  

  Note that moving row $i$ to a different address in memory requires an update
  to the ``row pointer'' in [[rows]];  this is why rows are stored along with
  their number.  When they need a new page, threads first try to grab an
  existing EMPTY page.  If there is none, a new page is allocated from the
  OS.  There are global doubly-linked lists of full and empty pages,
  protected by a lock, but these are infrequently accessed.
*/


int heap_config_get_PAGE_DATA_SIZE()
{
    return PAGE_DATA_SIZE;
}


/* provide an empty page */
static struct page_t *
heap_get_free_page(heapctx_t ctx)
{
        struct page_t *page = NULL;
        // #pragma omp critical(pagelist)
        omp_set_lock(&ctx->lock);
        {
                // try to grab it from the simply-linked list of empty pages.
                if (ctx->empty_pages != NULL) {
                        page = ctx->empty_pages->page;
                        ctx->empty_pages = ctx->empty_pages->next;
                        ctx->n_empty_pages--;
                } else {
                        ctx->n_pages++;   /* we will malloc() it, update count while still in critical section */
                }
        }
        omp_unset_lock(&ctx->lock);
        if (page == NULL) {
                // we must allocate a new page from the OS.
                page = malloc(sizeof(struct page_t));
                struct pagelist_t *item = malloc(sizeof(struct pagelist_t));
                page->list = item;
                item->page = page;
                page->i = ctx->n_pages;
        }
        page->ptr = 0;
        page->generation = ctx->current_generation;
        return page;
}

/* Provide the oldest full page with generation < current_generation, or NULL if none is available,
   and remove it from the doubly-linked list of full pages */
static struct page_t *
heap_get_full_page(heapctx_t ctx)
{
        struct page_t *page = NULL;
        //#pragma omp critical(pagelist)
        omp_set_lock(&ctx->lock);
        {
                struct pagelist_t *item = ctx->full_pages->next;
                if (item->page != NULL && item->page->generation < ctx->current_generation) {
                        page = item->page;
                        item->next->prev = item->prev;
                        item->prev->next = item->next;
                        ctx->n_full_pages--;
                }
        }
        omp_unset_lock(&ctx->lock);
        return page;
}


/* declare that the given page is empty */
static  void
heap_clear_page(heapctx_t ctx, struct page_t *page)
{
        struct pagelist_t *item = page->list;
        //#pragma omp critical(pagelist)
        omp_set_lock(&ctx->lock);
        {
                item->next = ctx->empty_pages;
                ctx->empty_pages = item;
                ctx->n_empty_pages++;
        }
        omp_unset_lock(&ctx->lock);
}

/* declare that the given page is full. Insert to the left of the list of full pages.
   The list is sorted (following next) by increasing generation. */
static void
heap_release_page(heapctx_t ctx, struct page_t *page)
{
        struct pagelist_t *list = page->list;
        omp_set_lock(&ctx->lock);
        {
                struct pagelist_t *target = ctx->full_pages->prev;
                while (target->page != NULL && page->generation < target->page->generation)
                        target = target->prev;
                list->next = target;
                list->prev = target->prev;
                list->next->prev = list;
                list->prev->next = list;
                ctx->n_full_pages++;
        }
        omp_unset_lock(&ctx->lock);
}

// set up the page linked lists
void
heap_setup(heapctx_t ctx)
{
        omp_init_lock(&ctx->lock);
        ctx->current_generation = 0;

        // setup the doubly-linked list of full pages.
        ctx->full_pages = &ctx->headnode;
        ctx->full_pages->page = NULL;
        ctx->full_pages->next = ctx->full_pages;
        ctx->full_pages->prev = ctx->full_pages;

        ctx->empty_pages = NULL;
        int T = omp_get_max_threads();
        ctx->active_page = malloc(T * sizeof(*ctx->active_page));
        ctx->heap_waste = malloc(T * sizeof(*ctx->heap_waste));

        #pragma omp parallel for
        for(int t = 0 ; t < T ; t++) {
            ctx->active_page[t] = heap_get_free_page(ctx);
            ctx->heap_waste[t] = 0;
        }
}

/* deallocate everything.  The memory is not released */
void
heap_reset(heapctx_t ctx)
{
        ctx->current_generation += 1;   // all full pages become "old"

        /* clear all full pages */
        struct page_t *page = heap_get_full_page(ctx);
        while (page != NULL) {
                heap_clear_page(ctx, page);
                page = heap_get_full_page(ctx);
        }

        /* clear active pages */
        int T = omp_get_max_threads ();
        for (int t = 0 ; t < T ; t++) {
                heap_clear_page(ctx, ctx->active_page[t]);
                ctx->active_page[t] = heap_get_free_page(ctx);
                ctx->heap_waste[t] = 0;
        }

        ctx->current_generation = 0;
}


/* release all memory */
void
heap_clear (heapctx_t ctx)
{
  /* clear active pages */
  int T = omp_get_max_threads ();
  for (int t = 0 ; t < T ; t++) {
    free(ctx->active_page[t]->list);
    free(ctx->active_page[t]);
  }

  /* clear empty pages */
  while (ctx->empty_pages != NULL) {
    struct pagelist_t *item = ctx->empty_pages;
    ctx->empty_pages = item->next;
    free(item->page);
    free(item);
  }

  /* clear full pages. 1. Locate dummy node */
  while (ctx->full_pages->page != NULL)
    ctx->full_pages = ctx->full_pages->next;

  // 2. Skip dummy node
  ctx->full_pages = ctx->full_pages->next;

  // 3. Walk list until dummy node is met again, free everything.
  while (ctx->full_pages->page != NULL) {
    struct pagelist_t *item = ctx->full_pages;
    ctx->full_pages = ctx->full_pages->next;
    free(item->page);
    free(item);
  }
  free(ctx->active_page);
  free(ctx->heap_waste);
}


/* Returns a pointer to allocated space holding a size-s array of typerow_t.
   Thread-safe. */
static inline typerow_t * heap_malloc (heapctx_t ctx, size_t s)
{
  ASSERT(s <= PAGE_DATA_SIZE);
  int t = omp_get_thread_num();
  struct page_t *page = ctx->active_page[t];
  // ASSERT(page != NULL);
  /* enough room in active page ?*/
  if (page->ptr + s >= PAGE_DATA_SIZE) {
        heap_release_page(ctx, page);
        page = heap_get_free_page(ctx);
        ctx->active_page[t] = page;
  }
  typerow_t *alloc = page->data + page->ptr;
  page->ptr += s;
  return alloc;
}


typerow_t *
heap_alloc_row (heapctx_t ctx, index_t i, size_t s)
{
  typerow_t *alloc = heap_malloc(ctx, s + 2);
  rowCell(alloc, 0) = i;
  rowCell(alloc, 1) = s;
  return alloc + 1;
}

void
heap_resize_last_row (heapctx_t ctx, typerow_t *row, index_t new_size)
{
  int t = omp_get_thread_num();
  struct page_t *page = ctx->active_page[t];
  index_t old_size = rowCell(row, 0);
  ASSERT(row + old_size + 1 == page->data + page->ptr);
  int delta = old_size - new_size;
  rowCell(row, 0) = new_size;
  page->ptr -= delta;
}


void
heap_destroy_row(heapctx_t ctx, typerow_t *row)
{
  int t = omp_get_thread_num();
  rowCell(row, -1) = (index_signed_t) -1;
  ctx->heap_waste[t] += rowCell(row, 0) + 2;
}

/* Copy non-garbage data to the active page of the current thread, then
   return the page to the freelist. Thread-safe.
   Warning: this may fill the current active page and release it. */
static int
collect_page(heapctx_t ctx, typerow_t **rows, struct page_t *page)
{
        int garbage = 0;
        int bot = 0;
        int top = page->ptr;
        typerow_t *data = page->data;
        while (bot < top) {
                index_signed_t i = rowCell(data, bot);
                typerow_t *old = data + bot + 1;
                index_t size = rowCell(old, 0);
                if (i == (index_signed_t) -1) {
                        garbage += size + 2;
                } else {
                        ASSERT(rows[i] == old);
                        typerow_t * new = heap_alloc_row(ctx, i, size);
                        memcpy(new, old, (size + 1) * sizeof(typerow_t));
                        setCell(new, -1, rowCell(old, -1), 0);
                        rows[i] = new;
                }
                bot += size + 2;
        }
        int t = omp_get_thread_num();
        ctx->heap_waste[t] -= garbage;
        heap_clear_page(ctx, page);
        return garbage;
}

static double
heap_waste_ratio(heapctx_t ctx)
{
        int T = omp_get_max_threads();
        long long total_waste = 0;
        for (int t = 0; t < T; t++)
                total_waste += ctx->heap_waste[t];
        double waste = ((double) total_waste) / (ctx->n_pages - ctx->n_empty_pages) / PAGE_DATA_SIZE;
        return waste;
}

void
heap_garbage_collection(heapctx_t ctx, typerow_t **rows)
{
        double waste = heap_waste_ratio(ctx);
        printf("Starting collection with %.0f%% of waste...", 100 * waste);
        fflush(stdout);

        // I don't want to collect pages just filled during the collection
        ctx->current_generation++;

        int i = 0;
        int initial_full_pages = ctx->n_full_pages;
        struct page_t *page;
        long long collected_garbage = 0;
        #pragma omp parallel reduction(+:i, collected_garbage) private(page)
        while ((page = heap_get_full_page(ctx)) != NULL) {
                collected_garbage += collect_page(ctx, rows, page);
                i++;
        }

        double page_ratio = (double) i / initial_full_pages;
        double recycling = 1 - heap_waste_ratio(ctx) / waste;
        if (i == 0)
                i = 1; // avoid division by zero
        printf("Examined %.0f%% of full pages, recycled %.0f%% of waste. %.0f%% of examined data was garbage\n",
        	100 * page_ratio, 100 * recycling, 100.0 * collected_garbage / i / PAGE_DATA_SIZE);
}
