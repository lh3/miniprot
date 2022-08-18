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
			r->vs = ((a[k]>>32) - mi->bo[ts]) << mi->opt.bbit;
			r->ve = ((a[k+n-1]>>32) - mi->bo[ts] + 1) << mi->opt.bbit;
			r->qs = (uint32_t)a[k];
			r->qe = (uint32_t)a[k+n-1];
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

void mp_set_parent(void *km, float mask_level, int mask_len, int n, mp_reg1_t *r, int sub_diff, int hard_mask_level) // and compute mp_reg1_t::subsc
{
	int i, j, k, *w;
	uint64_t *cov;
	if (n <= 0) return;
	for (i = 0; i < n; ++i) r[i].id = i;
	cov = (uint64_t*)kmalloc(km, n * sizeof(uint64_t));
	w = (int*)kmalloc(km, n * sizeof(int));
	w[0] = 0, r[0].parent = 0;
	for (i = 1, k = 1; i < n; ++i) {
		mp_reg1_t *ri = &r[i];
		int si = ri->qs, ei = ri->qe, n_cov = 0, uncov_len = 0;
		if (hard_mask_level) goto skip_uncov;
		for (j = 0; j < k; ++j) { // traverse existing primary hits to find overlapping hits
			mp_reg1_t *rp = &r[w[j]];
			int sj = rp->qs, ej = rp->qe;
			if (ej <= si || sj >= ei) continue;
			if (sj < si) sj = si;
			if (ej > ei) ej = ei;
			cov[n_cov++] = (uint64_t)sj<<32 | ej;
		}
		if (n_cov == 0) {
			goto set_parent_test; // no overlapping primary hits; then i is a new primary hit
		} else if (n_cov > 0) { // there are overlapping primary hits; find the length not covered by existing primary hits
			int j, x = si;
			radix_sort_mp64(cov, cov + n_cov);
			for (j = 0; j < n_cov; ++j) {
				if ((int)(cov[j]>>32) > x) uncov_len += (cov[j]>>32) - x;
				x = (int32_t)cov[j] > x? (int32_t)cov[j] : x;
			}
			if (ei > x) uncov_len += ei - x;
		}
skip_uncov:
		for (j = 0; j < k; ++j) { // traverse existing primary hits again
			mp_reg1_t *rp = &r[w[j]];
			int sj = rp->qs, ej = rp->qe, min, max, ol;
			if (ej <= si || sj >= ei) continue; // no overlap
			min = ej - sj < ei - si? ej - sj : ei - si;
			max = ej - sj > ei - si? ej - sj : ei - si;
			ol = si < sj? (ei < sj? 0 : ei < ej? ei - sj : ej - sj) : (ej < si? 0 : ej < ei? ej - si : ei - si); // overlap length; TODO: this can be simplified
			if ((float)ol / min - (float)uncov_len / max > mask_level && uncov_len <= mask_len) { // then this is a secondary hit
				int cnt_sub = 0, sci = ri->chn_sc;
				ri->parent = rp->parent;
				rp->subsc = rp->subsc > sci? rp->subsc : sci;
				if (ri->cnt >= rp->cnt) cnt_sub = 1;
				if (rp->p && ri->p && (rp->vid != ri->vid || rp->vs != ri->vs || rp->ve != ri->ve || ol != min)) { // the last condition excludes identical hits after DP
					sci = ri->p->dp_max;
					rp->p->dp_max2 = rp->p->dp_max2 > sci? rp->p->dp_max2 : sci;
					if (rp->p->dp_max - ri->p->dp_max <= sub_diff) cnt_sub = 1;
				}
				if (cnt_sub) ++rp->n_sub;
				break;
			}
		}
set_parent_test:
		if (j == k) w[k++] = i, ri->parent = i, ri->n_sub = 0;
	}
	kfree(km, cov);
	kfree(km, w);
}
