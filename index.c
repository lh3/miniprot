#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "mppriv.h"
#include "kalloc.h"
#include "kvec-km.h"
#include "kthread.h"

uint32_t *mp_idx_boff(const mp_ntdb_t *db, int32_t bbit, uint32_t *n_boff)
{
	int32_t i;
	int64_t boff = 0;
	uint32_t *bo;
	bo = Kmalloc(0, uint32_t, db->n_ctg * 2);
	for (i = 0; i < db->n_ctg; ++i) {
		bo[i<<1|0] = boff;
		boff += (db->ctg[i].len + (1<<bbit) - 1) >> bbit;
		bo[i<<1|1] = boff;
		boff += (db->ctg[i].len + (1<<bbit) - 1) >> bbit;
	}
	assert(boff < UINT32_MAX);
	*n_boff = boff;
	return bo;
}

void mp_idx_proc_orf(void *km, const uint8_t *seq, int64_t st, int64_t en, int32_t kmer, int32_t smer, int32_t bbit, int64_t boff, mp64_v *a)
{
	int64_t i;
	int32_t l;
	uint32_t x, mask_k = (1U<<kmer*4) - 1, mask_s = (1U<<smer*4) - 1;
	for (i = st, l = 0, x = 0; i < en; i += 3) {
		uint8_t codon = seq[i]<<4 | seq[i+1]<<2 | seq[i+2];
		x = (x<<4 | mp_tab_codon13[codon]) & mask_k;
		if (++l >= kmer) {
			int32_t sel = 0;
			uint32_t y = mp_hash32_mask(x, mask_k);
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
			if (sel)
				kv_push(uint64_t, km, *a, (uint64_t)y<<32 | ((st>>bbit) + boff));
		}
	}
}

void mp_idx_proc_seq(void *km, int64_t len, const uint8_t *seq, int32_t min_aa_len, int32_t kmer, int32_t smer, int32_t bbit, int64_t boff, mp64_v *a)
{
	uint8_t codon[3], p;
	int64_t i, j, e[3], k[3], l[3];
	kv_resize(uint64_t, km, *a, len>>2);
	a->n = 0;
	for (p = 0; p < 3; ++p)
		e[p] = -1, k[p] = l[p] = 0, codon[p] = 0;
	for (i = 0, p = 0; i < len; ++i, ++p) {
		if (p == 3) p = 0;
		if (seq[i] < 4) {
			codon[p] = (codon[p] << 2 | seq[i]) & 0x3f;
			if (++l[p] >= 3) {
				uint8_t aa = mp_tab_codon[(uint8_t)codon[p]];
				if (aa >= 20) {
					if (k[p] >= min_aa_len)
						mp_idx_proc_orf(km, seq, e[p] + 1 - k[p] * 3, e[p] + 1, kmer, smer, bbit, boff, a);
					k[p] = l[p] = 0, e[p] = -1;
				} else e[p] = i, ++k[p];
			}
		} else {
			if (k[p] >= min_aa_len)
				mp_idx_proc_orf(km, seq, e[p] + 1 - k[p] * 3, e[p] + 1, kmer, smer, bbit, boff, a);
			k[p] = l[p] = 0, e[p] = -1;
		}
	}
	for (p = 0; p < 3; ++p)
		if (k[p] >= min_aa_len)
			mp_idx_proc_orf(km, seq, e[p] + 1 - k[p] * 3, e[p] + 1, kmer, smer, bbit, boff, a);
	radix_sort_mp64(a->a, a->a + a->n);
	for (i = 1, j = 0; i < a->n; ++i)
		if (a->a[j] != a->a[i])
			a->a[++j] = a->a[i];
	a->n = ++j;
}

typedef struct {
	const mp_idxopt_t *io;
	const mp_idx_t *mi;
	void **km;
	mp64_v *a;
} worker_aux_t;

static void build_worker(void *data, long j, int tid)
{
	worker_aux_t *aux = (worker_aux_t*)data;
	const mp_ntdb_t *nt = aux->mi->nt;
	const mp_idxopt_t *io = aux->io;
	void *km = aux->km[tid];
	int64_t len;
	uint8_t *seq;
	seq = Kmalloc(km, uint8_t, nt->ctg[j>>1].len);
	len = mp_ntseq_get(nt, j>>1, 0, -1, j&1, seq);
	mp_idx_proc_seq(km, len, seq, io->min_aa_len, io->kmer, io->smer, io->bbit, aux->mi->bo[j], &aux->a[j]);
	kfree(km, seq);
}

static void build_bidx(const mp_idxopt_t *io, mp_idx_t *mi, const mp64_v *a)
{
	int32_t i, j, n_a = mi->nt->n_ctg * 2;
	int64_t tmp;
	uint32_t n_kmer = 1U << io->kmer*4;
	mi->ki = Kcalloc(0, int64_t, n_kmer);
	for (i = 0; i < n_a; ++i) { // count
		const mp64_v *b = &a[i];
		for (j = 0; j < b->n; ++j)
			++mi->ki[b->a[j]>>32];
	}
	for (i = 0, mi->n_kb = 0; i < n_kmer; ++i)
		mi->n_kb += mi->ki[i], mi->ki[i] = mi->n_kb - mi->ki[i];
	mi->kb = Kmalloc(0, uint32_t, mi->n_kb);
	for (i = 0; i < n_a; ++i) {
		const mp64_v *b = &a[i];
		for (j = 0; j < b->n; ++j)
			mi->kb[mi->ki[b->a[j]>>32]++] = (uint32_t)b->a[j];
	}
	for (i = 0, tmp = 0; i < n_kmer; ++i) {
		int64_t t = tmp;
		tmp = mi->ki[i], mi->ki[i] -= tmp - t;
	}
}

mp_idx_t *mp_idx_build(const mp_idxopt_t *io, const char *fn, int32_t n_threads)
{
	mp_ntdb_t *nt;
	mp_idx_t *mi;
	worker_aux_t aux;
	int32_t i;

	nt = mp_ntseq_read(fn);
	if (nt == 0) return 0;

	mi = Kcalloc(0, mp_idx_t, 1);
	mi->bbit = io->bbit, mi->kmer = io->kmer, mi->smer = io->smer;
	mi->nt = nt;
	mi->bo = mp_idx_boff(mi->nt, io->bbit, &mi->n_block);

	memset(&aux, 0, sizeof(aux));
	aux.io = io;
	aux.mi = mi;
	aux.km = Kcalloc(0, void*, n_threads);
	aux.a = Kcalloc(0, mp64_v, nt->n_ctg * 2);
	for (i = 0; i < n_threads; ++i)
		aux.km[i] = km_init();
	kt_for(n_threads, build_worker, &aux, nt->n_ctg * 2);
	build_bidx(io, mi, aux.a);
	fprintf(stderr, "[M::%s] %d blocks and %ld kmer-block pairs\n", __func__, mi->n_block, (long)mi->n_kb);

	for (i = 0; i < n_threads; ++i)
		km_destroy(aux.km[i]);
	free(aux.a); free(aux.km);
	return mi;
}

void mp_idx_destroy(mp_idx_t *mi)
{
	if (mi == 0) return;
	mp_ntseq_destroy(mi->nt);
	free(mi->ki); free(mi->bo); free(mi->kb);
	free(mi);
}

void mp_idx_dump(FILE *fp, const mp_idx_t *mi)
{
	int32_t x[3];
	x[0] = mi->bbit, x[1] = mi->kmer, x[2] = mi->smer;
	fwrite(MP_MAGIC, 1, 4, fp);
	fwrite(x, 4, 3, fp);
	fwrite(&mi->n_kb, 8, 1, fp);
	mp_ntseq_dump(fp, mi->nt);
	fwrite(mi->ki, 8, 1U<<mi->kmer*4, fp);
	fwrite(mi->kb, 4, mi->n_kb, fp);
}

mp_idx_t *mp_idx_restore(FILE *fp)
{
	char magic[4];
	int32_t x[3];
	mp_idx_t *mi;
	fread(magic, 1, 4, fp);
	if (strncmp(magic, MP_MAGIC, 4) != 0)
		return 0;
	mi = Kcalloc(0, mp_idx_t, 1);
	fread(x, 4, 3, fp);
	mi->bbit = x[0], mi->kmer = x[1], mi->smer = x[2];
	fread(&mi->n_kb, 8, 1, fp);
	mi->nt = mp_ntseq_restore(fp);
	mi->ki = Kmalloc(0, int64_t, 1U<<mi->kmer*4);
	mi->kb = Kmalloc(0, uint32_t, mi->n_kb);
	fread(mi->ki, 8, 1U<<mi->kmer*4, fp);
	fread(mi->kb, 4, mi->n_kb, fp);
	return mi;
}
