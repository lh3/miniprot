#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include "mppriv.h"
#include "kthread.h"
#include "kvec-km.h"
#include "kalloc.h"
#include "kseq.h"

struct mp_tbuf_s {
	void *km;
	int rep_len;
};

mp_tbuf_t *mp_tbuf_init(void)
{
	mp_tbuf_t *b;
	b = Kcalloc(0, mp_tbuf_t, 1);
	if (!(mp_dbg_flag & MP_DBG_NO_KALLOC))
		b->km = km_init();
	return b;
}

void mp_tbuf_destroy(mp_tbuf_t *b)
{
	if (b == 0) return;
	km_destroy(b->km);
	free(b);
}

static void mp_refine_reg(void *km, const mp_idx_t *mi, const mp_mapopt_t *opt, const char *aa, int32_t l_aa, mp_reg1_t *r, int32_t extl, int32_t extr)
{
	const mp_idxopt_t *io = &mi->opt;
	int32_t i, j, k, n_u, max_sc = 0, max_i, kmer = opt->kmer2, is_splice = !(opt->flag&MP_F_NO_SPLICE);
	int64_t as, ae, l_nt, n_a, ctg_len = mi->nt->ctg[r->vid>>1].len;
	uint8_t *nt;
	uint64_t *a, *u;
	mp64_v sd = {0,0,0}, sd_aa = {0,0,0};

	as = r->vs > extl? r->vs - extl : 0;
	ae = r->ve + extr < ctg_len? r->ve + extr : ctg_len;
	nt = Kmalloc(km, uint8_t, ae - as);
	l_nt = mp_ntseq_get(mi->nt, r->vid>>1, r->vid&1? ctg_len - ae : as, r->vid&1? ctg_len - as : ae, r->vid&1, nt);
	mp_sketch_nt4(km, nt, l_nt, io->min_aa_len, kmer, 0, 0, 0, &sd);
	kfree(km, nt);
	mp_sketch_prot(km, aa, l_aa, kmer, 0, &sd_aa);
	for (i = 0; i < sd_aa.n; ++i)
		kv_push(uint64_t, km, sd, sd_aa.a[i] | 1ULL<<31);
	kfree(km, sd_aa.a);
	radix_sort_mp64(sd.a, sd.a + sd.n);

	for (k = 0, i = 1, n_a = 0; i <= sd.n; ++i) { // precalculate n_a
		if (i == sd.n || sd.a[k]>>32 != sd.a[i]>>32) {
			int32_t n1, n2;
			for (j = k; j < i; ++j)
				if (sd.a[j]>>31&1) break;
			n1 = j - k, n2 = i - k - n1;
			if (n1 > 0 && n2 > 0 && n1 * n2 <= opt->max_ava)
				n_a += n1 * n2;
			k = i;
		}
	}
	a = Kmalloc(km, uint64_t, n_a);
	for (k = 0, i = 1, n_a = 0; i <= sd.n; ++i) { // populate a[]
		if (i == sd.n || sd.a[k]>>32 != sd.a[i]>>32) {
			int32_t n1, n2;
			for (j = k; j < i; ++j)
				if (sd.a[j]>>31&1) break;
			n1 = j - k, n2 = i - k - n1;
			if (n1 > 0 && n2 > 0 && n1 * n2 <= opt->max_ava) {
				int32_t i1, i2;
				for (i1 = k; i1 < k + n1; ++i1)
					for (i2 = k + n1; i2 < i; ++i2)
						a[n_a++] = (uint64_t)((uint32_t)sd.a[i1])<<32 | ((uint32_t)sd.a[i2]<<1>>1);
			}
			k = i;
		}
	}
	kfree(km, sd.a);

	radix_sort_mp64(a, a + n_a);
	a = mp_chain(opt->max_intron, opt->max_gap, opt->bw, opt->max_chn_max_skip, opt->max_chn_iter, opt->min_chn_cnt, opt->min_chn_sc, opt->chn_coef_log, is_splice, kmer, 0, n_a, a, &n_u, &u, km);
	if (n_u == 0) {
		r->cnt = 0, r->off = -1, r->a = 0;
		return;
	}

	max_sc = u[0]>>32, max_i = 0;
	for (i = 1; i < n_u; ++i)
		if (max_sc < u[i]>>32)
			max_sc = u[i]>>32, max_i = i;
	for (i = k = 0; i < n_u; ++i) {
		if (i == max_i) break;
		k += (uint32_t)u[i];
	}

	n_a = (int32_t)u[max_i];
	memmove(a, a + k, sizeof(*a) * n_a);
	r->a = a = krelocate(km, a, sizeof(*a) * n_a);
	r->chn_sc = u[max_i]>>32;
	r->cnt = n_a, r->off = 0;
	r->qs = (uint32_t)a[0] - (kmer - 1);
	r->qe = (uint32_t)a[n_a-1] + 1;
	r->vs = as + (a[0]>>32) + 1 - 3 * kmer;
	r->ve = as + (a[n_a-1]>>32) + 1;
	for (i = 0; i < n_a; ++i)
		a[i] = ((a[i]>>32) + as - r->vs) << 32 | a[i]<<32>>32;
	r->chn_sc_ungap = mp_cal_chn_sc_ungap(r->cnt, r->a, kmer);
	kfree(km, u);
}

static void mp_dbg_chain(const mp_idx_t *mi, int32_t n_reg, const mp_reg1_t *reg, const uint64_t *a, int32_t bbit, const char *label)
{
	int32_t i, k;
	for (i = 0; i < n_reg; ++i) {
		const mp_reg1_t *r = &reg[i];
		for (k = 0; k < r->cnt; ++k) {
			uint64_t ak = a[r->off + k];
			uint64_t off = bbit > 0? (long)((ak>>32) - mi->bo[r->vid]) << bbit : r->vs + (ak>>32);
			fprintf(stderr, "%s\t%d\t%ld\t%s\t%c\t%ld\t%d\n", label, i, (long)(ak>>32), mi->nt->ctg[r->vid>>1].name, "+-"[r->vid&1], (long)off, (uint32_t)ak);
		}
	}
}

static int32_t mp_cal_max_occ(void *km, const mp_idx_t *mi, int32_t n, const uint64_t *a)
{
	int32_t i, n_kmer;
	uint64_t *cnt, q25, q75;
	cnt = Kmalloc(km, uint64_t, n);
	n_kmer = mp_n_bucket(&mi->opt);
	for (i = 0; i < n; ++i) {
		int64_t en = (a[i]>>32) + 1 < n_kmer? mi->ki[(a[i]>>32) + 1] : mi->n_kb;
		cnt[i] = en - mi->ki[a[i]>>32];
	}
	radix_sort_mp64(cnt, cnt + n);
	q25 = cnt[(int64_t)(n * .25 + .499)];
	q75 = cnt[(int64_t)(n * .75 + .499)];
	kfree(km, cnt);
	return (int32_t)(q75 + (q75 - q25) * 1.5 + 10.); // threshold the outliers, the boxplot way (except for +10)
}

mp_reg1_t *mp_map(const mp_idx_t *mi, int qlen, const char *seq, int *n_reg, mp_tbuf_t *b, const mp_mapopt_t *opt, const char *qname)
{
	void *km = b->km;
	const mp_idxopt_t *io = &mi->opt;
	mp64_v sd = {0,0,0};
	mp_reg1_t *reg;
	int32_t i, n_u, max_occ, is_splice = !(opt->flag&MP_F_NO_SPLICE);
	int64_t k, n_a = 0, n_kmer;
	uint64_t *a, *u;

	*n_reg = 0;
	n_kmer = mp_n_bucket(&mi->opt);
	mp_sketch_prot(km, seq, qlen, io->kmer, io->mod_bit, &sd);
	radix_sort_mp64(sd.a, sd.a + sd.n);

	if (sd.n >= 8) {
		max_occ = mp_cal_max_occ(km, mi, sd.n, sd.a);
		if (max_occ > opt->max_occ) max_occ = opt->max_occ;
	} else max_occ = opt->max_occ;

	for (i = 0; i < sd.n; ++i) {
		int64_t en = (sd.a[i]>>32) + 1 < n_kmer? mi->ki[(sd.a[i]>>32) + 1] : mi->n_kb;
		int64_t n = en - mi->ki[sd.a[i]>>32];
		if (n <= max_occ) n_a += n;
	}
	a = Kmalloc(km, uint64_t, n_a);
	for (i = 0, k = 0; i < sd.n; ++i) {
		int64_t j, st = mi->ki[sd.a[i]>>32];
		int64_t en = (sd.a[i]>>32) + 1 < n_kmer? mi->ki[(sd.a[i]>>32) + 1] : mi->n_kb;
		if (en - st <= max_occ)
			for (j = st; j < en; ++j)
				a[k++] = (uint64_t)mi->kb[j] << 32 | (uint32_t)sd.a[i];
	}
	kfree(km, sd.a);
	radix_sort_mp64(a, a + n_a);

	if (mp_dbg_flag & MP_DBG_ANCHOR) {
		for (k = 0; k < n_a; ++k) {
			i = mp_idx_block2pos(mi, a[k]>>32);
			fprintf(stderr, "X\t%ld\t%s\t%c\t%ld\t%d\n", (long)(a[k]>>32), mi->nt->ctg[i>>1].name, "+-"[i&1], (long)((a[k]>>32) - mi->bo[i]) << mi->opt.bbit, (uint32_t)a[k]);
		}
	}

	if (!(opt->flag & MP_F_NO_PRE_CHAIN) && is_splice) {
		int32_t w = 1<<mi->opt.bbit;
		a = mp_chain(w, w, w, opt->max_chn_max_skip, opt->max_chn_iter, 2, 0, opt->chn_coef_log, is_splice, mi->opt.kmer, mi->opt.bbit, n_a, a, &n_u, &u, km);
		for (i = 0, n_a = 0; i < n_u; ++i) n_a += (uint32_t)u[i];
		kfree(km, u);
		radix_sort_mp64(a, a + n_a);
	}

	a = mp_chain(opt->max_intron, opt->max_gap, opt->bw, opt->max_chn_max_skip, opt->max_chn_iter, opt->min_chn_cnt, opt->min_chn_sc, opt->chn_coef_log,
				 is_splice, mi->opt.kmer, mi->opt.bbit, n_a, a, &n_u, &u, km);
	reg = mp_reg_gen_from_block(0, mi, n_u, u, a, n_reg);
	kfree(km, u);
	mp_sort_reg(km, n_reg, reg);
	mp_set_parent(km, opt->mask_level, opt->mask_len, *n_reg, reg, mi->opt.kmer, 0);
	mp_select_sub(km, opt->pri_ratio * opt->pri_ratio, mi->opt.kmer * 2, opt->best_n, n_reg, reg);

	if (mp_dbg_flag & MP_DBG_CHAIN)
		mp_dbg_chain(mi, *n_reg, reg, a, mi->opt.bbit, "Y1");

	if (!(mp_dbg_flag & MP_DBG_NO_REFINE)) {
		int32_t nr = 0;
		uint64_t *ext;
		ext = mp_cal_max_ext(km, *n_reg, reg, a, 100, opt->max_ext);
		for (i = 0; i < *n_reg; ++i) {
			mp_refine_reg(km, mi, opt, seq, qlen, &reg[i], ext[i]>>32, (int32_t)ext[i]);
			if (reg[i].cnt > 0)
				reg[nr++] = reg[i];
		}
		*n_reg = nr;
		kfree(km, ext);
		kfree(km, a);
		a = mp_collate_a(km, *n_reg, reg);
		mp_sort_reg(km, n_reg, reg);
		//mp_dbg_chain(mi, *n_reg, reg, a, 0, "Y2");
		mp_set_parent(km, opt->mask_level, opt->mask_len, *n_reg, reg, mi->opt.kmer, 0);
		mp_select_sub(km, opt->pri_ratio * opt->pri_ratio, mi->opt.kmer * 2, opt->best_n, n_reg, reg);
	}
	if (!(opt->flag & MP_F_NO_ALIGN)) {
		int32_t k;
		for (i = k = 0; i < *n_reg; ++i) {
			mp_align(km, opt, mi, qlen, seq, &reg[i]);
			if (reg[i].p) reg[k++] = reg[i];
		}
		*n_reg = k;
		mp_sort_reg(km, n_reg, reg);
		mp_select_multi_exon(*n_reg, reg, opt->io);
		mp_set_parent(km, opt->mask_level, opt->mask_len, *n_reg, reg, mi->opt.kmer, 0);
		mp_select_sub(km, opt->pri_ratio, mi->opt.kmer * 2, opt->best_n, n_reg, reg);
	}
	kfree(km, a);
	return reg;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int32_t n_threads;
	int64_t id;
	const mp_mapopt_t *opt;
	mp_bseq_file_t *fp;
	const mp_idx_t *mi;
	kstring_t str;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int32_t n_seq;
	mp_bseq1_t *seq;
	int32_t *n_reg;
	mp_reg1_t **reg;
	mp_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *s = (step_t*)_data;
	mp_bseq1_t *seq = &s->seq[i];
	if (mp_dbg_flag & MP_DBG_QNAME)
		fprintf(stderr, "QR\t%s\t%d\t%d\n", s->seq[i].name, s->seq[i].l_seq, tid);
	s->reg[i] = mp_map(s->p->mi, seq->l_seq, seq->seq, &s->n_reg[i], s->buf[tid], s->p->opt, seq->name);
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int32_t i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
        s = Kcalloc(0, step_t, 1);
		s->seq = mp_bseq_read(p->fp, p->opt->mini_batch_size, 0, &s->n_seq);
		if (s->seq) {
			s->p = p;
			s->buf = Kcalloc(0, mp_tbuf_t*, p->n_threads);
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = mp_tbuf_init();
			s->n_reg = Kcalloc(0, int32_t, s->n_seq);
			s->reg = Kcalloc(0, mp_reg1_t*, s->n_seq);
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: map
		kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_seq);
		return in;
    } else if (step == 2) { // step 2: output
        step_t *s = (step_t*)in;
		for (i = 0; i < p->n_threads; ++i) mp_tbuf_destroy(s->buf[i]);
		free(s->buf);
		for (i = 0; i < s->n_seq; ++i) {
			int32_t best_sc = -1, j, n_out;
			if (s->n_reg[i] > 0)
				best_sc = s->reg[i][0].p? s->reg[i][0].p->dp_max : s->reg[i][0].chn_sc;
			for (j = 0, n_out = 0; j < s->n_reg[i] && j < p->opt->out_n; ++j) {
				const mp_reg1_t *r = &s->reg[i][j];
				int32_t sc = r->p? r->p->dp_max : r->chn_sc;
				if (sc <= 0 || sc < (double)best_sc * p->opt->out_sim) continue;
				if (r->qe - r->qs < (double)s->seq[i].l_seq * p->opt->out_cov) continue;
				mp_write_output(&p->str, 0, p->mi, &s->seq[i], r, p->opt, ++p->id, j + 1);
				fwrite(p->str.s, 1, p->str.l, stdout);
				++n_out;
			}
			if (n_out == 0) {
				mp_write_output(&p->str, 0, p->mi, &s->seq[i], 0, p->opt, 0, 0);
				fwrite(p->str.s, 1, p->str.l, stdout);
			}
			for (j = 0; j < s->n_reg[i]; ++j) {
				free(s->reg[i][j].feat);
				free(s->reg[i][j].p);
			}
			free(s->reg[i]);
			free(s->seq[i].seq); free(s->seq[i].name);
			if (s->seq[i].comment) free(s->seq[i].comment);
		}
		free(s->reg); free(s->n_reg); free(s->seq);
		if (mp_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, mp_realtime(), mp_cputime() / mp_realtime(), s->n_seq);
		free(s);
	}
    return 0;
}

int32_t mp_map_file(const mp_idx_t *idx, const char *fn, const mp_mapopt_t *opt, int n_threads)
{
	pipeline_t pl;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.fp = mp_bseq_open(fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.mi = idx;
	pl.n_threads = n_threads > 1? n_threads : 1;
	if (opt->flag & MP_F_GFF) puts("##gff-version 3");
	kt_pipeline(2, worker_pipeline, &pl, 3);
	free(pl.str.s);
	mp_bseq_close(pl.fp);
	return 0;
}
