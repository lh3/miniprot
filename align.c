#include <stdio.h>
#include <assert.h>
#include "mppriv.h"
#include "nasw.h"

static void mp_filter_seed(int32_t cnt, uint64_t *a, int32_t max_aa_dist, int32_t min_cnt)
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
			for (; i < j; ++i)
				a[i] |= 1ULL<<31;
			--i;
		}
	}
}

typedef struct {
	int32_t n, m;
	uint32_t *c;
} mp_cigar_t;

static void mp_align_seq(void *km, const mp_mapopt_t *opt, int32_t nlen, const uint8_t *nseq, int32_t alen, const char *aseq, mp_cigar_t *cigar)
{
	if (nlen == alen * 3 && alen <= opt->kmer2) {
		cigar->c = ns_push_cigar(km, &cigar->n, &cigar->m, cigar->c, NS_CIGAR_M, alen);
	} else {
		int32_t i;
		ns_opt_t ns_opt;
		ns_rst_t rst;
		ns_opt_init(&ns_opt);
		ns_opt.flag |= NS_F_CIGAR;
		ns_opt.go = opt->go, ns_opt.ge = opt->ge, ns_opt.io = opt->io, ns_opt.fs = opt->fs, ns_opt.nc = opt->nc, ns_opt.sc = opt->mat;
		ns_rst_init(&rst);
		ns_global_gs16(km, (const char*)nseq, nlen, aseq, alen, &ns_opt, &rst);
		printf("%d\t%d\t%d\t%d\n", nlen, alen, rst.score, rst.n_cigar);
		for (i = 0; i < rst.n_cigar; ++i)
			cigar->c = ns_push_cigar(km, &cigar->n, &cigar->m, cigar->c, rst.cigar[i]&0xf, rst.cigar[i]>>4);
	}
}

static void mp_extra_cal(mp_reg1_t *r, const int8_t *mat)
{
	int32_t k, nl = 0, al = 0;
	mp_extra_t *e = r->p;
	for (k = 0; k < e->n_cigar; ++k) {
		int32_t op = e->cigar[k]&0xf, len = e->cigar[k]>>4;
		if (op == NS_CIGAR_M) {
			nl += len * 3, al += len;
		} else if (op == NS_CIGAR_I) {
			al += len;
		} else if (op == NS_CIGAR_D) {
			nl += len * 3;
		} else if (op == NS_CIGAR_F) {
			nl += len;
		} else if (op == NS_CIGAR_G) {
			nl += len, ++al;
		} else if (op == NS_CIGAR_N) {
			nl += len;
		} else if (op == NS_CIGAR_U) {
			nl += len, ++al;
		} else if (op == NS_CIGAR_V) {
			nl += len, ++al;
		}
	}
	printf("%d == %lld\t%d == %d\n", nl, r->ve - r->vs, al, r->qe - r->qs);
}

static void mp_extra_gen(void *km, mp_reg1_t *r, mp_cigar_t *cigar, int32_t score)
{
	int32_t cap;
	kfree(km, r->a);
	r->a = 0;
	cap = cigar->n + sizeof(mp_extra_t) / 4;
	r->p = Kcalloc(0, mp_extra_t, cap); // allocate globally, not from km
	r->p->dp_max = r->p->dp_max2 = score;
	r->p->n_cigar = r->p->m_cigar = cigar->n;
	memcpy(r->p->cigar, cigar->c, cigar->n * sizeof(uint32_t));
	kfree(km, cigar->c);
}

void mp_align(void *km, const mp_mapopt_t *opt, const mp_idx_t *mi, int32_t len, const char *aa, mp_reg1_t *r)
{
	int32_t i, i0, ne0, ae0;
	int64_t as, ae, ctg_len, l_nt;
	uint8_t *nt;
	mp_cigar_t cigar = {0,0,0};

	assert(r->cnt > 0);
	mp_filter_seed(r->cnt, r->a, 3, 3);

	ctg_len = mi->nt->ctg[r->vid>>1].len;
	as = r->vs > opt->max_ext? r->vs - opt->max_ext : 0;
	ae = r->ve + opt->max_ext < ctg_len? r->ve + opt->max_ext : ctg_len;
	nt = Kmalloc(km, uint8_t, ae - as);
	l_nt = mp_ntseq_get(mi->nt, r->vid>>1, r->vid&1? ctg_len - ae : as, r->vid&1? ctg_len - as : ae, r->vid&1, nt);
	assert(l_nt == ae - as);

	for (i = 0; i < r->cnt; ++i)
		if (r->a[i]>>31&1) break;
	i0 = i;
	#if 1
	if (i0 < r->cnt) {
		cigar.c = ns_push_cigar(km, &cigar.n, &cigar.m, cigar.c, NS_CIGAR_M, opt->kmer2);
		ne0 = (r->a[i0]>>32) + 1, ae0 = ((int32_t)r->a[i0]<<1>>1) + 1;
	}
	for (i = i0 + 1; i < r->cnt; ++i) {
		int32_t ne1, ae1;
		if (!(r->a[i]>>31&1)) continue;
		ne1 = (r->a[i]>>32) + 1, ae1 = ((int32_t)r->a[i]<<1>>1) + 1;
		mp_align_seq(km, opt, ne1 - ne0, &nt[ne0 + r->vs - as], ae1 - ae0, &aa[ae0], &cigar);
		i0 = i, ne0 = ne1, ae0 = ae1;
	}
	#else
	mp_align_seq(km, opt, r->ve - r->vs, &nt[r->vs - as], r->qe - r->qs, aa, &cigar);
	#endif
	kfree(km, nt);
	printf("%lld\n", as);
	for (i = 0; i < cigar.n; ++i) printf("%d%c", cigar.c[i]>>4, NS_CIGAR_STR[cigar.c[i]&0xf]); putchar('\n');
	mp_extra_gen(km, r, &cigar, 0);
	mp_extra_cal(r, opt->mat);
}
