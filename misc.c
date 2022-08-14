#include "mppriv.h"
#include "ksort.h"

#define sort_key_64(x) (x)
KRADIX_SORT_INIT(mp64, uint64_t, sort_key_64, 8)

int32_t mp_verbose = 3;

void mp_start(void)
{
	mp_make_tables(MP_CODON_STD);
	mp_realtime();
}
