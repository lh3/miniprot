#include <stdlib.h>
#include "mppriv.h"
#include "kalloc.h"

mp_reg1_t *mp_reg_gen_from_block(void *km, const mp_idx_t *mi, int32_t n_u, const uint64_t *u, const uint64_t *a, int32_t *n_reg)
{
	int32_t i, k, nr;
	mp_reg1_t *reg;
	reg = Kcalloc(km, mp_reg1_t, n_u);
	for (i = k = nr = 0; i < n_u; ++i) {
		uint32_t n = (uint32_t)u[i];
		int32_t ts, te;
		mp_reg1_t *r = &reg[nr++];
		r->off = k, r->cnt = n;
		ts = mp_idx_block2pos(mi, a[k]>>32);
		te = mp_idx_block2pos(mi, a[k+n-1]>>32);
		if (ts == te) { // on the same contig
			r->vid = ts;
			r->st = ((a[k]>>32) - mi->bo[ts]) << mi->opt.bbit;
			r->en = ((a[k+n-1]>>32) - mi->bo[ts] + 1) << mi->opt.bbit;
		} else { // on different contigs
			fprintf(stderr, "ERROR: not implemented yet\n");
			--nr;
		}
		r->chn_sc = u[i] >> 32;
		k += n;
	}
	*n_reg = nr;
	return reg;
}
