#ifndef MPPRIV_H
#define MPPRIV_H

#include "miniprot.h"

#ifndef kroundup64
#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#endif

#ifdef __cplusplus
extern "C" {
#endif

static inline uint32_t mp_hash32_mask(uint32_t key, uint32_t mask)
{
	key += ~(key << 15) & mask;
	key ^=  (key >> 10);
	key +=  (key << 3) & mask;
	key ^=  (key >> 6);
	key += ~(key << 11) & mask;
	key ^=  (key >> 16);
	return key;
}

void radix_sort_mp64(uint64_t *st, uint64_t *en);

double mp_realtime(void);
double mp_cputime(void);
long mp_peakrss(void);
char *mp_strdup(const char *src);

void mp_idx_proc_seq(void *km, int64_t len, const uint8_t *seq, int32_t min_aa_len, int32_t kmer, int32_t smer, int32_t bbit, int64_t boff, mp64_v *a);

#ifdef __cplusplus
}
#endif

#endif
