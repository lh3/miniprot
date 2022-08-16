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
	if (!(mp_dbg_flag & 1)) b->km = km_init();
	return b;
}

void mp_tbuf_destroy(mp_tbuf_t *b)
{
	if (b == 0) return;
	km_destroy(b->km);
	free(b);
}

mp_reg1_t *mp_map(const mp_idx_t *mi, int qlen, const char *seq, int *n_reg, mp_tbuf_t *b, const mp_mapopt_t *opt, const char *qname)
{
	void *km = b->km;
	const mp_idxopt_t *io = &mi->opt;
	mp64_v sd = {0,0,0};
	mp_reg1_t *reg;

	*n_reg = 0;
	mp_sketch_prot(km, seq, qlen, io->kmer, io->smer, &sd);

	int32_t i;
	int64_t k, n_a = 0;
	for (i = 0; i < sd.n; ++i) { // TODO: sorting might help to reduce cache misses, but probably doesn't matter in practice
		int64_t n = mi->ki[(sd.a[i]>>32) + 1] - mi->ki[sd.a[i]>>32];
		if (n <= opt->max_occ) n_a += n;
	}
	uint64_t *a;
	a = Kmalloc(km, uint64_t, n_a);
	for (i = 0, k = 0; i < sd.n; ++i) {
		int64_t j, st = mi->ki[sd.a[i]>>32], en = mi->ki[(sd.a[i]>>32) + 1];
		if (en - st <= opt->max_occ)
			for (j = st; j < en; ++j)
				a[k++] = (uint64_t)mi->kb[j] << 32 | (uint32_t)sd.a[i];
	}
	radix_sort_mp64(a, a + n_a);

	int32_t n_u;
	uint64_t *u;
	a = mp_chain(opt->max_intron, opt->max_gap, opt->bw, 25, 1000000, opt->min_chn_cnt, 0, 1, mi->opt.kmer, mi->opt.bbit, n_a, a, &n_u, &u, km);
	reg = mp_reg_gen_from_block(0, mi, n_u, u, a, n_reg);
	/*
	for (i = 0; i < *n_reg; ++i) {
		mp_reg1_t *r = &reg[i];
		fprintf(stderr, "%s\t%ld\t%ld\t%s\t%d\t%c\n", mi->nt->ctg[r->vid>>1].name, (long)r->vs, (long)r->ve, qname, r->chn_sc, "+-"[r->vid&1]);
	}
	if (0) {
		fprintf(stderr, "NC\t%d\n", n_u);
		for (i = 0, n_a = 0; i < n_u; ++i) {
			n_a += (int32_t)u[i];
			fprintf(stderr, "CN\t%d\t%d\t%d\n", i, (int32_t)(u[i]>>32), (int32_t)u[i]);
		}
	}
	*/
	//for (k = 0; k < n_a; ++k) printf("%s\t%d\t%d\t%d\n", qname, (uint32_t)(a[k]>>32)<<mi->opt.bbit, (uint32_t)(a[k]>>32), (uint32_t)a[k]);
	kfree(km, a);
	kfree(km, sd.a);
	return reg;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int32_t n_threads;
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
		int32_t j;
		for (i = 0; i < p->n_threads; ++i) mp_tbuf_destroy(s->buf[i]);
		free(s->buf);
		for (i = 0; i < s->n_seq; ++i) {
			for (j = 0; j < s->n_reg[i]; ++j) {
				p->str.l = 0;
				mp_write_paf(&p->str, p->mi, &s->seq[i], &s->reg[i][j]);
				fwrite(p->str.s, 1, p->str.l, stdout);
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
	kt_pipeline(2, worker_pipeline, &pl, 3);
	free(pl.str.s);
	mp_bseq_close(pl.fp);
	return 0;
}
