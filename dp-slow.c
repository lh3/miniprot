#include "mppriv.h"
#include "kalloc.h"

typedef struct {
	void *km;
	int32_t len, asize;
	uint8_t *seq;
	int8_t *qp;
	uint8_t mem[];
} mp_dps_qp_t;

mp_dps_qp_t *mp_dps_qprof(void *km, const char *seq, int32_t len, int32_t asize, const int8_t *mat)
{
	int32_t i, j, k;
	mp_dps_qp_t *q;
	q = (mp_dps_qp_t*)kmalloc(km, sizeof(mp_dps_qp_t) * (len + len * asize));
	q->len = len;
	q->asize = asize;
	q->seq = q->mem;
	q->qp = (int8_t*)(q->seq + len);
	for (j = 0; j < len; ++j)
		q->seq[j] = mp_tab_aa20[(uint8_t)seq[j]];
	for (i = 0, k = 0; i < asize; ++i) {
		const int8_t *p = &mat[asize * i];
		for (j = 0; j < len; ++j)
			q->qp[k++] = p[q->seq[j]];
	}
	return q;
}
