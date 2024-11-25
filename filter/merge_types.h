#ifndef MERGE_TYPES_H_
#define MERGE_TYPES_H_

#include "typedefs.h"

#ifndef FOR_DL
#define typerow_t index_t
#define cmp_typerow_t cmp_index
#else
#define typerow_t ideal_merge_t
#define cmp_typerow_t cmp_ideal_merge
#endif

#define col_weight_t unsigned char

#endif