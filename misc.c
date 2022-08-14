#include "ksort.h"

#define sort_key_64(x) (x)
KRADIX_SORT_INIT(mp64, uint64_t, sort_key_64, 8)
