#include "mppriv.h"
#include "ksort.h"

#define sort_key_64(x) (x)
KRADIX_SORT_INIT(mp64, uint64_t, sort_key_64, 8)

#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(mp128x, mp128_t, sort_key_128x, 8)

int32_t mp_verbose = 3, mp_dbg_flag = 0;

void mp_start(void)
{
	ns_make_tables(MP_CODON_STD);
	mp_realtime();
}

char *mp_strdup(const char *src)
{
	size_t len;
	char *dst;
	len = strlen(src);
	dst = (char*)malloc(len + 1);
	memcpy(dst, src, len + 1);
	return dst;
}
