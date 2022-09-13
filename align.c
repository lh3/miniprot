#include <stdio.h>
#include <assert.h>
#include "mppriv.h"
#include "nasw.h"

static void mp_filter_seed(int32_t cnt, uint64_t *a, int32_t max_aa_dist, int32_t min_cnt, int32_t kmer2, int32_t trim_back)
{
	int32_t i, j;
	for (i = 0; i < cnt; ++i) {
		for (j = i + 1; j < cnt; ++j) {
			int32_t x0, y0, x1, y1;
			x0 = a[j-1]>>32, y0 = (int32_t)a[j-1];
			x1 = a[j]>>32,   y1 = (int32_t)a[j];
			if ((x1 - x0) % 3 != 0 || x1 - x0 > max_aa_dist * 3 || y1 - y0 > max_aa_dist)
				break;
		}
		if (j - i >= min_cnt) {
			int32_t k, t = (int32_t)a[j-1];
			for (k = j - 2; k >= i; --k)
				if (t - (int32_t)a[k] >= trim_back)
					break;
			t = (int32_t)a[i] + 1 - kmer2;
			for (; i < k; ++i)
				if ((int32_t)a[i] + 1 - t >= trim_back)
					break;
			for (; i <= k; ++i)
				a[i] |= 1ULL<<31;
			i = j - 1;
		}
	}
}

static inline int32_t mp_score_ungapped(int32_t alen, const uint8_t *nseq, const char *aseq, int32_t asize, const int8_t *mat)
{
	int32_t i, j, score = 0;
	for (i = 0, j = 0; i < alen; i += 3, ++j) {
		uint8_t nt_aa, aa_aa, codon = nseq[i]<<4 | nseq[i+1]<<2 | nseq[i+2];
		nt_aa = nseq[i] > 3 || nseq[i+1] > 3 || nseq[i+2] > 3? ns_tab_aa20['X'] : ns_tab_codon[codon];
		aa_aa = ns_tab_aa20[(uint8_t)aseq[j]];
		score += mat[nt_aa * asize + aa_aa];
	}
	return score;
}

typedef struct {
	int32_t n, m;
	uint32_t *c;
} mp_cigar_t;

static inline void mp_map2ns_opt(const mp_mapopt_t *mo, ns_opt_t *no)
{
	ns_opt_init(no);
	no->go = mo->go, no->ge = mo->ge, no->io = mo->io, no->fs = mo->fs, no->nc = mo->nc, no->sc = mo->mat;
	no->end_bonus = mo->end_bonus;
}

static int32_t mp_align_seq(void *km, const mp_mapopt_t *opt, int32_t nlen, const uint8_t *nseq, int32_t alen, const char *aseq, mp_cigar_t *cigar)
{
	int32_t i;
	if (nlen == alen * 3 && alen <= opt->kmer2) {
		cigar->c = ns_push_cigar(km, &cigar->n, &cigar->m, cigar->c, NS_CIGAR_M, alen);
		return mp_score_ungapped(alen, nseq, aseq, opt->asize, opt->mat);
	} else {
		ns_opt_t ns_opt;
		ns_rst_t rst;
		mp_map2ns_opt(opt, &ns_opt);
		ns_opt.flag |= NS_F_CIGAR;
		ns_rst_init(&rst);
		ns_global_gs16(km, (const char*)nseq, nlen, aseq, alen, &ns_opt, &rst);
		//printf("%d\t%d\t%d\t%d\n", nlen, alen, rst.score, rst.n_cigar);
		for (i = 0; i < rst.n_cigar; ++i)
			cigar->c = ns_push_cigar(km, &cigar->n, &cigar->m, cigar->c, rst.cigar[i]&0xf, rst.cigar[i]>>4);
		kfree(km, rst.cigar);
		return rst.score;
	}
}

static void mp_extra_cal(mp_reg1_t *r, const mp_mapopt_t *opt, const uint8_t *nt, const char *aa)
{
	int32_t k, i, j, l, nl = 0, al = 0;
	mp_extra_t *e = r->p;
	e->clen = e->n_iden = e->n_plus = e->dp_max = 0, e->dist_stop = -1;
	for (k = 0; k < e->n_cigar; ++k) {
		int32_t op = e->cigar[k]&0xf, len = e->cigar[k]>>4, len3 = len * 3;
		if (op == NS_CIGAR_M) {
			for (i = nl, j = al, l = 0; l < len; ++l, ++j, i += 3) {
				uint8_t nt_aa, aa_aa, codon = nt[i]<<4 | nt[i+1]<<2 | nt[i+2];
				int32_t s;
				nt_aa = nt[i] > 3 || nt[i+1] > 3 || nt[i+2] > 3? ns_tab_aa20['X'] : ns_tab_codon[codon];
				aa_aa = ns_tab_aa20[(uint8_t)aa[j]];
				s = opt->mat[nt_aa * opt->asize + aa_aa];
				e->n_iden += (nt_aa == aa_aa);
				e->n_plus += (s > 0);
				e->dp_max += s;
			}
			nl += len3, al += len, e->clen += len3;
		} else if (op == NS_CIGAR_I) {
			e->dp_max -= opt->go + opt->ge * len;
			al += len, e->clen += len3;
		} else if (op == NS_CIGAR_D) {
			e->dp_max -= opt->go + opt->ge * len;
			nl += len3, e->clen += len3;
		} else if (op == NS_CIGAR_F) {
			e->dp_max -= opt->fs;
			nl += len, e->clen += len;
		} else if (op == NS_CIGAR_G) {
			e->dp_max -= opt->fs;
			nl += len, ++al, e->clen += 3;
		} else if (op == NS_CIGAR_N) {
			nl += len;
		} else if (op == NS_CIGAR_U || op == NS_CIGAR_V) {
			uint8_t n1, n2, n3, codon, nt_aa, aa_aa;
			int32_t s;
			if (op == NS_CIGAR_U) n1 = nt[nl], n2 = nt[nl + len - 2], n3 = nt[nl + len - 1];
			else n1 = nt[nl], n2 = nt[nl + 1], n3 = nt[nl + len - 1];
			codon = n1<<4 | n2<<2 | n3;
			nt_aa = n1 > 3 || n2 > 3 || n3 > 3? ns_tab_aa20['X'] : ns_tab_codon[codon];
			aa_aa = ns_tab_aa20[(uint8_t)aa[al]];
			s = opt->mat[nt_aa * opt->asize + aa_aa];
			//printf("%c:%c, score=%d\n", ns_tab_aa_i2c[nt_aa], ns_tab_aa_i2c[aa_aa], s);
			e->n_iden += (nt_aa == aa_aa);
			e->n_plus += (s > 0);
			e->dp_max += s;
			nl += len, ++al, e->clen += 3;
		}
	}
	if (nl != r->ve - r->vs || al != r->qe - r->qs) {
		fprintf(stderr, "BUG! %d == %d? %d == %d? ", nl, (int)(r->ve - r->vs), al, r->qe - r->qs);
		for (k = 0; k < e->n_cigar; ++k)
			fprintf(stderr, "%d%c", e->cigar[k]>>4, NS_CIGAR_STR[e->cigar[k]&0xf]);
		fputc('\n', stderr);
	}
	assert(nl == r->ve - r->vs);
	assert(al == r->qe - r->qs);
	//printf("aa_score=%d, glen=%d, clen=%d, n_iden=%d, n_plus=%d\n", e->aa_score, e->glen, e->clen, e->n_iden, e->n_plus);
}

static void mp_extra_gen(void *km, mp_reg1_t *r, mp_cigar_t *cigar, int32_t score)
{
	int32_t cap;
	cap = cigar->n + sizeof(mp_extra_t) / 4;
	r->p = Kcalloc(0, mp_extra_t, cap); // allocate globally, not from km
	r->p->dp_score = score;
	r->p->n_cigar = r->p->m_cigar = cigar->n;
	memcpy(r->p->cigar, cigar->c, cigar->n * sizeof(uint32_t));
	kfree(km, cigar->c);
}

static int32_t mp_extra_stop(const mp_reg1_t *r, const uint8_t *nt, int64_t as, int64_t ae)
{
	int64_t j;
	for (j = r->ve; j + 2 < ae; j += 3) {
		int32_t i = j - as;
		uint8_t codon = nt[i]<<4 | nt[i+1]<<2 | nt[i+2];
		uint8_t aa = nt[i] > 3 || nt[i+1] > 3 || nt[i+2] > 3? ns_tab_aa20['X'] : ns_tab_codon[codon];
		if (aa == 20) return j - r->ve;
	}
	return -1;
}

static int32_t mp_extra_start(const mp_reg1_t *r, const uint8_t *nt, int64_t as, int64_t ae)
{
	int64_t j;
	for (j = r->vs; j >= as && j + 2 < ae; j -= 3) {
		int32_t i = j - as;
		uint8_t codon = nt[i]<<4 | nt[i+1]<<2 | nt[i+2];
		uint8_t aa = nt[i] > 3 || nt[i+1] > 3 || nt[i+2] > 3? ns_tab_aa20['X'] : ns_tab_codon[codon];
		if (aa == 20) break;
		if (aa == 12) return r->vs - j;
	}
	return -1;
}

void mp_align(void *km, const mp_mapopt_t *opt, const mp_idx_t *mi, int32_t len, const char *aa, mp_reg1_t *r)
{
	int32_t i, i0, ne0 = 0, ae0 = 0, score = 0, extl, extr;
	int64_t as, ae, ctg_len, vs0, l_nt;
	uint8_t *nt;
	mp_cigar_t cigar = {0,0,0};

	assert(r->cnt > 0);
	mp_filter_seed(r->cnt, r->a, 3, 3, opt->kmer2, opt->kmer2 + 1);
	for (i = 0; i < r->cnt; ++i)
		if (r->a[i]>>31&1) break;
	if (i == r->cnt) { // all filtered; FIXME: we need to filter it later; not implemented yet
		r->cnt = 0;
		return;
	}
	i0 = i;

	extl = extr = opt->max_ext;
	if (r->qs >= 10) extl = opt->max_intron/2;
	if (len - r->qe >= 10) extr = opt->max_intron/2;
	ctg_len = mi->nt->ctg[r->vid>>1].len;
	as = r->vs > extl? r->vs - extl : 0;
	ae = r->ve + extr < ctg_len? r->ve + extr : ctg_len;
	nt = Kmalloc(km, uint8_t, ae - as);
	l_nt = mp_ntseq_get_by_v(mi->nt, r->vid, as, ae, nt);
	assert(l_nt == ae - as);
	vs0 = r->vs;
	//fprintf(stderr, "X\t%s\t%c\t%ld\t%ld\n", mi->nt->ctg[r->vid>>1].name, "+-"[r->vid&1], (long)(r->vid&1? ctg_len - ae : as), (long)(r->vid&1? ctg_len - as : ae));

	{ // left extension
		int64_t vs1 = vs0 + (r->a[i0]>>32) + 1;
		int32_t as1 = ((int32_t)r->a[i0]<<1>>1) + 1;
		ns_opt_t ns_opt;
		ns_rst_t rst;
		mp_map2ns_opt(opt, &ns_opt);
		ns_opt.flag |= NS_F_EXT_LEFT;
		ns_rst_init(&rst);
		ns_global_gs16(km, (const char*)nt, vs1 - as, aa, as1, &ns_opt, &rst);
		r->vs = vs1 - rst.nt_len;
		r->qs = as1 - rst.aa_len;
		ne0 = r->vs - vs0;
		ae0 = r->qs;
	}

	#if 1
	for (i = i0; i < r->cnt; ++i) {
		int32_t ne1, ae1;
		if (!(r->a[i]>>31&1)) continue;
		ne1 = (r->a[i]>>32) + 1, ae1 = ((int32_t)r->a[i]<<1>>1) + 1;
		score += mp_align_seq(km, opt, ne1 - ne0, &nt[ne0 + vs0 - as], ae1 - ae0, &aa[ae0], &cigar);
		i0 = i, ne0 = ne1, ae0 = ae1;
	}
	r->ve = ne0 + vs0, r->qe = ae0;
	#else
	score = mp_align_seq(km, opt, r->ve - r->vs, &nt[r->vs - as], r->qe - r->qs, aa, &cigar);
	#endif

	if (r->qe < len && r->ve < ae) { // right extension
		ns_opt_t ns_opt;
		ns_rst_t rst;
		mp_map2ns_opt(opt, &ns_opt);
		ns_opt.flag |= NS_F_EXT_RIGHT;
		ns_rst_init(&rst);
		ns_global_gs16(km, (const char*)&nt[r->ve - as], ae - r->ve, &aa[r->qe], len - r->qe, &ns_opt, &rst);
		score += mp_align_seq(km, opt, rst.nt_len, &nt[r->ve - as], rst.aa_len, &aa[r->qe], &cigar);
		r->ve += rst.nt_len, r->qe += rst.aa_len;
	}

	//for (i = 0; i < cigar.n; ++i) printf("%d%c", cigar.c[i]>>4, NS_CIGAR_STR[cigar.c[i]&0xf]); putchar('\n');
	mp_extra_gen(km, r, &cigar, score);
	mp_extra_cal(r, opt, &nt[r->vs - as], &aa[r->qs]);
	r->p->dist_stop  = mp_extra_stop(r, nt, as, ae);
	r->p->dist_start = mp_extra_start(r, nt, as, ae);
	kfree(km, nt);
}
