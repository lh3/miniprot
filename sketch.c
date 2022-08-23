#include <assert.h>
#include <stdio.h>
#include "mppriv.h"
#include "kalloc.h"
#include "kvec-km.h"

static inline uint32_t mp_hash32_mask(uint32_t key, uint32_t mask)
{
	key  = (key + ~(key << 15)) & mask;
	key ^=  (key >> 10);
	key  =  (key + (key <<  3)) & mask;
	key ^=  (key >> 6);
	key  = (key + ~(key << 11)) & mask;
	key ^=  (key >> 16);
	return key;
}

void mp_sketch_prot(void *km, const char *seq, int32_t len, int32_t kmer, int32_t smer, mp64_v *a)
{
	int32_t i, j, l = 0;
	uint32_t x = 0, y, mask_k = (1U<<kmer*4) - 1, mask_s = (1U<<smer*4) - 1;
	a->n = 0;
	for (i = 0; i < len; ++i) {
		int32_t c = ns_tab_aa13[(uint8_t)seq[i]];
		if (c < 14) {
			x = (x<<4 | c) & mask_k;
			if (++l >= kmer) {
				int32_t sel = 0;
				y = mp_hash32_mask(x, mask_k);
				if (kmer > smer) {
					uint32_t m = UINT32_MAX;
					for (j = 0; j < (kmer - smer) << 2; j += 2) {
						uint32_t z = y>>j & mask_s;
						m = m < z? m : z;
					}
					sel = (m == (y>>(kmer-smer)*2 & mask_s));
				} else sel = 1;
				if (sel) kv_push(uint64_t, km, *a, (uint64_t)y<<32 | i);
			}
		} else x = 0, l = 0;
	}
}

void mp_sketch_clean_orf(void *km, const uint8_t *seq, int64_t st, int64_t en, int32_t kmer, int32_t smer, int32_t bbit, int64_t boff, mp64_v *a)
{
	int64_t i;
	int32_t l;
	uint32_t x, mask_k = (1U<<kmer*4) - 1, mask_s = (1U<<smer*4) - 1;
	for (i = st, l = 0, x = 0; i < en; i += 3) {
		uint8_t codon = seq[i]<<4 | seq[i+1]<<2 | seq[i+2];
		x = (x<<4 | ns_tab_codon13[codon]) & mask_k;
		if (++l >= kmer) {
			int32_t sel = 0;
			uint32_t y;
			y = mp_hash32_mask(x, mask_k);
			assert(y <= mask_k);
			if (kmer > smer) {
				int32_t j;
				uint32_t m = UINT32_MAX;
				for (j = 0; j < (kmer - smer) << 2; j += 2) {
					uint32_t z = y>>j & mask_s;
					m = m < z? m : z;
				}
				sel = (m == (y>>(kmer-smer)*2 & mask_s));
			} else sel = 1;
			if (sel) kv_push(uint64_t, km, *a, (uint64_t)y<<32 | (((i+2)>>bbit) + boff));
		}
	}
}

void mp_sketch_nt4(void *km, const uint8_t *seq, int64_t len, int32_t min_aa_len, int32_t kmer, int32_t smer, int32_t bbit, int64_t boff, mp64_v *a)
{
	uint8_t codon, p;
	int64_t i, j, l, e[3], k[3];
	kv_resize(uint64_t, km, *a, len>>2);
	a->n = 0;
	for (p = 0; p < 3; ++p)
		e[p] = -1, k[p] = 0;
	for (i = 0, p = 1, codon = 0, l = 0; i < len; ++i, ++p) {
		if (p == 3) p = 0;
		if (seq[i] < 4) {
			codon = (codon << 2 | seq[i]) & 0x3f;
			if (++l >= 3) {
				uint8_t aa = ns_tab_codon[codon];
				if (aa >= 20) {
					if (k[p] >= min_aa_len)
						mp_sketch_clean_orf(km, seq, e[p] + 1 - k[p] * 3, e[p] + 1, kmer, smer, bbit, boff, a);
					k[p] = 0, e[p] = -1;
				} else e[p] = i, ++k[p];
			}
		} else {
			if (k[p] >= min_aa_len)
				mp_sketch_clean_orf(km, seq, e[p] + 1 - k[p] * 3, e[p] + 1, kmer, smer, bbit, boff, a);
			k[p] = 0, e[p] = -1, l = 0, codon = 0;
		}
	}
	for (p = 0; p < 3; ++p)
		if (k[p] >= min_aa_len)
			mp_sketch_clean_orf(km, seq, e[p] + 1 - k[p] * 3, e[p] + 1, kmer, smer, bbit, boff, a);
	radix_sort_mp64(a->a, a->a + a->n);
	for (i = 1, j = 0; i < a->n; ++i)
		if (a->a[j] != a->a[i])
			a->a[++j] = a->a[i];
	a->n = ++j;
}
