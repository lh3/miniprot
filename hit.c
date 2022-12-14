#include <stdlib.h>
#include <assert.h>
#include "mppriv.h"
#include "kalloc.h"

static int32_t mp_cal_chn_sc_ungap_approx(const mp_reg1_t *r, const uint64_t *a, int32_t kmer)
{
	int32_t i, x = kmer;
	for (i = 1; i < r->cnt; ++i) {
		const uint64_t a0 = a[r->off + i - 1], a1 = a[r->off + i];
		int32_t dq = (int32_t)a1 - (int32_t)a0;
		x += dq < kmer? dq : kmer;
		if (a1>>32 == a0>>32) x += MP_BLOCK_BONUS;
	}
	return x;
}

int32_t mp_cal_chn_sc_ungap(int32_t n_a, const uint64_t *a, int32_t kmer)
{
	int32_t i, x = kmer;
	for (i = 1; i < n_a; ++i) {
		const uint64_t a0 = a[i - 1], a1 = a[i];
		int32_t dq = (int32_t)a1 - (int32_t)a0, dr3 = (a1>>32) - (a0>>32);
		int32_t dr = dr3 / 3, q = dr3 - dr * 3, dg;
		dg = dq < dr? dq : dr;
		if (dq >= dr && q != 0) --x;
		else x += dg < kmer? dg : kmer;
	}
	return x;
}

mp_reg1_t *mp_reg_gen_from_block(void *km, const mp_idx_t *mi, int32_t n_u, const uint64_t *u, const uint64_t *a, int32_t *n_reg)
{
	int32_t i, k, nr;
	mp_reg1_t *reg;
	reg = Kcalloc(km, mp_reg1_t, n_u);
	for (i = k = nr = 0; i < n_u; ++i) {
		uint32_t n = (uint32_t)u[i];
		int32_t ts, te, is, ie;
		mp_reg1_t *r = &reg[nr++];
		r->off = k, r->cnt = n;
		is = k, ie = k + n - 1;
		ts = mp_idx_block2pos(mi, a[is]>>32);
		te = mp_idx_block2pos(mi, a[ie]>>32);
		assert(ts <= te);
		if (ts == te) { // on the same contig
			r->vid = ts;
		} else { // on different contigs
			int32_t j, js, je;
			for (j = k; j < k + n; ++j)
				if (a[j]>>32 >= mi->bo[ts+1])
					break;
			assert(j < k + n);
			js = j;
			for (j = k + n - 1; j >= js; --j)
				if (a[j]>>32 < mi->bo[te])
					break;
			je = j + 1;
			if (js - k > k + n - je) {
				r->vid = ts, ie = js - 1;
			} else {
				r->vid = te, is = je;
			}
			//fprintf(stderr, "SPLIT: %d<%d %d<=%d<=%d<=%d\n", ts, te, k, js, je, k + n);
		}
		r->vs = ((a[is]>>32) - mi->bo[r->vid]) << mi->opt.bbit;
		r->ve = ((a[ie]>>32) - mi->bo[r->vid] + 1) << mi->opt.bbit;
		r->qs = (uint32_t)a[is];
		r->qe = (uint32_t)a[ie] + 0;
		r->chn_sc = ts == te? u[i]>>32 : (uint32_t)((double)(u[i]>>32) * (ie - is + 1) / n + .499);
		r->chn_sc_ungap = mp_cal_chn_sc_ungap_approx(r, a, mi->opt.kmer);
		k += n;
	}
	*n_reg = nr;
	return reg;
}

uint64_t *mp_collate_a(void *km, int32_t n_reg, mp_reg1_t *reg)
{
	uint64_t *a;
	int64_t n_a;
	int32_t i;
	for (i = 0, n_a = 0; i < n_reg; ++i)
		n_a += reg[i].cnt;
	a = Kmalloc(km, uint64_t, n_a);
	for (i = 0, n_a = 0; i < n_reg; ++i) {
		mp_reg1_t *r = &reg[i];
		r->off = n_a;
		memcpy(&a[n_a], r->a, r->cnt * sizeof(uint64_t));
		kfree(km, r->a);
		r->a = &a[n_a];
		n_a += r->cnt;
	}
	return a;
}

void mp_sort_reg(void *km, int *n_regs, mp_reg1_t *r)
{
	int32_t i, n_aux, n = *n_regs, has_cigar = 0, no_cigar = 0;
	mp128_t *aux;
	mp_reg1_t *t;

	if (n <= 1) return;
	aux = (mp128_t*)kmalloc(km, n * 16);
	t = (mp_reg1_t*)kmalloc(km, n * sizeof(mp_reg1_t));
	for (i = n_aux = 0; i < n; ++i) {
		if (r[i].cnt > 0) { // squeeze out elements with cnt==0 (soft deleted)
			int score;
			if (r[i].p) score = r[i].p->dp_max, has_cigar = 1;
			else score = r[i].chn_sc, no_cigar = 1;
			aux[n_aux].x = (uint64_t)score << 32 | r[i].hash;
			aux[n_aux++].y = i;
		} else if (r[i].p) {
			free(r[i].p);
			r[i].p = 0;
		}
	}
	assert(has_cigar + no_cigar == 1);
	radix_sort_mp128x(aux, aux + n_aux);
	for (i = n_aux - 1; i >= 0; --i)
		t[n_aux - 1 - i] = r[aux[i].y];
	memcpy(r, t, sizeof(mp_reg1_t) * n_aux);
	*n_regs = n_aux;
	kfree(km, aux);
	kfree(km, t);
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

void mp_sync_regs(void *km, int n_regs, mp_reg1_t *regs) // keep mp_reg1_t::{id,parent} in sync; also reset id
{
	int *tmp, i, max_id = -1, n_tmp;
	if (n_regs <= 0) return;
	for (i = 0; i < n_regs; ++i) // NB: doesn't work if mp_reg1_t::id is negative
		max_id = max_id > regs[i].id? max_id : regs[i].id;
	n_tmp = max_id + 1;
	tmp = (int*)kmalloc(km, n_tmp * sizeof(int));
	for (i = 0; i < n_tmp; ++i) tmp[i] = -1;
	for (i = 0; i < n_regs; ++i)
		if (regs[i].id >= 0) tmp[regs[i].id] = i;
	for (i = 0; i < n_regs; ++i) {
		mp_reg1_t *r = &regs[i];
		r->id = i;
		if (r->parent == MP_PARENT_TMP_PRI)
			r->parent = i;
		else if (r->parent >= 0 && tmp[r->parent] >= 0)
			r->parent = tmp[r->parent];
		else r->parent = MP_PARENT_UNSET;
	}
	kfree(km, tmp);
}

void mp_select_sub(void *km, float pri_ratio, int min_diff, int best_n, int *n_, mp_reg1_t *r)
{
	if (pri_ratio > 0.0f && *n_ > 0) {
		int i, k, n = *n_, n_2nd = 0;
		for (i = k = 0; i < n; ++i) {
			int p = r[i].parent;
			int sci = r[i].p? r[i].p->dp_max : r[i].chn_sc;
			int scp = r[p].p? r[p].p->dp_max : r[p].chn_sc;
			if (p == i) { // primary or inversion
				r[k++] = r[i];
			} else if ((sci >= scp * pri_ratio || sci + min_diff >= scp) && n_2nd < best_n) {
				if (!(r[i].qs == r[p].qs && r[i].qe == r[p].qe && r[i].vid == r[p].vid && r[i].vs == r[p].vs && r[i].ve == r[p].ve)) // not identical hits
					r[k++] = r[i], ++n_2nd;
				else if (r[i].p) free(r[i].p);
			} else if (r[i].p) free(r[i].p);
		}
		if (k != n) mp_sync_regs(km, k, r); // removing hits requires sync()
		*n_ = k;
	}
}

void mp_select_multi_exon(int32_t n, mp_reg1_t *r, int32_t single_penalty)
{
	int32_t i;
	mp_reg1_t t;
	if (n < 2 || r[0].n_exon != 1) return;
	for (i = 1; i < n; ++i)
		if (r[i].n_exon >= 2)
			break;
	if (i == n) return;
	if (r[0].p == 0 || r[i].p == 0) return;
	if (r[0].p->dp_max < r[i].p->dp_max + single_penalty)
		t = r[0], r[0] = r[i], r[i] = t;
}

uint64_t *mp_cal_max_ext(void *km, int32_t n_reg, mp_reg1_t *reg, const uint64_t *a, int32_t min_ext, int32_t max_ext)
{
	uint64_t *ext, *b;
	int32_t i;
	if (n_reg <= 0) return 0;
	b = Kmalloc(km, uint64_t, n_reg);
	for (i = 0; i < n_reg; ++i) {
		mp_reg1_t *r = &reg[i];
		b[i] = (a[r->off]>>32)<<32 | i;
	}
	radix_sort_mp64(b, b + n_reg);
	ext = Kmalloc(km, uint64_t, n_reg);
	for (i = 0; i < n_reg; ++i) {
		int32_t left = max_ext, right = max_ext, j = (int32_t)b[i];
		mp_reg1_t *r = &reg[j];
		if (i > 0) {
			mp_reg1_t *q = &reg[(int32_t)b[i-1]];
			if (q->vid == r->vid && q->qe >= r->qs) {
				left = r->vs - q->ve < max_ext? r->vs - q->ve : max_ext;
				left = left > min_ext? left : min_ext;
			}
		}
		if (i < n_reg - 1) {
			mp_reg1_t *q = &reg[(int32_t)b[i+1]];
			if (q->vid == r->vid && r->qe >= q->qs) {
				right = q->vs - r->ve < max_ext? q->vs - r->ve : max_ext;
				right = right > min_ext? right : min_ext;
			}
		}
		ext[j] = (uint64_t)left<<32 | right;
	}
	kfree(km, b);
	return ext;
}
