#ifndef MST_H_
#define MST_H_

#include "filter_config.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

int minimalSpanningTree(int *start, int *end, int m, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX]);

#ifdef __cplusplus
}
#endif

#endif	/* MST_H_ */
