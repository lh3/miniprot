#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <fcntl.h>
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
	bo = Kmalloc(0, uint32_t, db->n_ctg * 2 + 1);
	for (i = 0; i < db->n_ctg; ++i) {
		bo[i<<1|0] = boff;
		boff += (db->ctg[i].len + (1<<bbit) - 1) >> bbit;
		bo[i<<1|1] = boff;
		boff += (db->ctg[i].len + (1<<bbit) - 1) >> bbit;
	}
	assert(boff < UINT32_MAX);
	bo[db->n_ctg*2] = *n_boff = boff;
	return bo;
}

static int32_t mp_idx_block2pos_core(const uint32_t *a, uint32_t n_a, uint32_t b)
{
	int32_t st = 0, en = n_a - 1;
	if (b >= a[n_a]) return -1;
	while (st <= en) {
		int32_t mid = (st + en) / 2;
		if (a[mid] <= b && b < a[mid+1]) return mid;
		else if (b < a[mid]) en = mid - 1;
		else st = mid + 1;
	}
	return -2;
}

int32_t mp_idx_block2pos(const mp_idx_t *mi, uint32_t b)
{
	return mp_idx_block2pos_core(mi->bo, mi->nt->n_ctg * 2, b);
}

typedef struct {
	const mp_idx_t *mi;
	void **km, **km_tmp;
	mp64_v *a;
} worker_aux_t;

static void build_worker(void *data, long j, int tid)
{
	worker_aux_t *aux = (worker_aux_t*)data;
	const mp_ntdb_t *nt = aux->mi->nt;
	const mp_idxopt_t *io = &aux->mi->opt;
	void *km = aux->km[tid], *km_tmp = aux->km_tmp[tid];
	int64_t len;
	uint8_t *seq;
	mp64_v a = {0,0,0};
	seq = Kmalloc(km_tmp, uint8_t, nt->ctg[j>>1].len);
	len = mp_ntseq_get(nt, j>>1, 0, -1, j&1, seq);
	mp_sketch_nt4(km_tmp, seq, len, io->min_aa_len, io->kmer, io->mod_bit, io->bbit, aux->mi->bo[j], &a);
	aux->a[j].m = aux->a[j].n = a.n;
	aux->a[j].a = Kmalloc(km, uint64_t, a.n);
	memcpy(aux->a[j].a, a.a, a.n * sizeof(uint64_t));
	kfree(km_tmp, a.a);
	kfree(km_tmp, seq);
}

static void build_bidx(const mp_idxopt_t *io, mp_idx_t *mi, const mp64_v *a)
{
	int32_t i, j, n_a = mi->nt->n_ctg * 2;
	int64_t tmp;
	uint32_t n_kmer;
	n_kmer = mp_n_bucket(io);
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

mp_idx_t *mp_idx_build(const char *fn, const mp_idxopt_t *io, int32_t n_threads)
{
	mp_ntdb_t *nt;
	mp_idx_t *mi;
	worker_aux_t aux;
	int32_t i;

	nt = mp_ntseq_read(fn);
	if (nt == 0) return 0;

	mi = Kcalloc(0, mp_idx_t, 1);
	mi->opt = *io;
	mi->nt = nt;
	mi->bo = mp_idx_boff(mi->nt, io->bbit, &mi->n_block);
	if (mp_verbose >= 3)
		fprintf(stderr, "[M::%s@%.3f*%.2f] %d blocks\n", __func__, mp_realtime(), mp_percent_cpu(), mi->n_block);

	memset(&aux, 0, sizeof(aux));
	aux.mi = mi;
	aux.km = Kcalloc(0, void*, n_threads);
	aux.km_tmp = Kcalloc(0, void*, n_threads);
	aux.a = Kcalloc(0, mp64_v, nt->n_ctg * 2);
	for (i = 0; i < n_threads; ++i) {
		aux.km[i] = km_init();
		aux.km_tmp[i] = km_init();
	}
	kt_for(n_threads, build_worker, &aux, nt->n_ctg * 2);
	if (mp_verbose >= 3)
		fprintf(stderr, "[M::%s@%.3f*%.2f] collected syncmers\n", __func__, mp_realtime(), mp_percent_cpu());
	build_bidx(io, mi, aux.a);
	if (mp_verbose >= 3)
		fprintf(stderr, "[M::%s@%.3f*%.2f] %ld kmer-block pairs\n", __func__, mp_realtime(), mp_percent_cpu(), (long)mi->n_kb);

	for (i = 0; i < n_threads; ++i) {
		km_destroy(aux.km[i]);
		km_destroy(aux.km_tmp[i]);
	}
	free(aux.a); free(aux.km);
	return mi;
}

void mp_idx_print_stat(const mp_idx_t *mi, int32_t max_occ)
{
	int64_t tot = 0, sum_large = 0;
	uint32_t i, n, n_occupied = 0, n_large = 0;
	n = mp_n_bucket(&mi->opt);
	for (i = 0; i < n - 1; ++i) {
		int64_t c = mi->ki[i+1] - mi->ki[i];
		if (c > 0) ++n_occupied;
		if (c > max_occ) ++n_large, sum_large += c;
		else tot += c;
	}
	fprintf(stderr, "[M::%s] %d distinct k-mers; mean occ of infrequent k-mers: %.2f; %d frequent k-mers accounting for %ld occurrences\n", __func__,
			n_occupied, (double)tot / (n_occupied - n_large), n_large, (long)sum_large);
}

void mp_idx_destroy(mp_idx_t *mi)
{
	if (mi == 0) return;
	mp_ntseq_destroy(mi->nt);
	free(mi->ki); free(mi->bo); free(mi->kb);
	free(mi);
}

/*
 * Index I/O
 */

static int64_t mp_idx_is_idx(const char *fn)
{
	int fd, is_idx = 0;
	int64_t ret, off_end;
	char magic[4];

	if (strcmp(fn, "-") == 0) return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0) return -1; // error
#ifdef WIN32
	if ((off_end = _lseeki64(fd, 0, SEEK_END)) >= 4) {
		_lseeki64(fd, 0, SEEK_SET);
#else
	if ((off_end = lseek(fd, 0, SEEK_END)) >= 4) {
		lseek(fd, 0, SEEK_SET);
#endif // WIN32
		ret = read(fd, magic, 4);
		if (ret == 4 && strncmp(magic, MP_IDX_MAGIC, 4) == 0)
			is_idx = 1;
	}
	close(fd);
	return is_idx? off_end : 0;
}

int mp_idx_dump(const char *fn, const mp_idx_t *mi)
{
	FILE *fp;
	fp = strcmp(fn, "-") == 0? stdout : fopen(fn, "wb");
	if (fp == 0) return -1;
	fwrite(MP_IDX_MAGIC, 1, 4, fp);
	fwrite(&mi->opt, sizeof(mi->opt), 1, fp);
	fwrite(&mi->n_kb, 8, 1, fp);
	mp_ntseq_dump(fp, mi->nt);
	fwrite(mi->ki, 8, mp_n_bucket(&mi->opt), fp);
	fwrite(mi->kb, 4, mi->n_kb, fp);
	if (fp != stdout) fclose(fp);
	return 0;
}

mp_idx_t *mp_idx_restore(const char *fn)
{
	FILE *fp;
	char magic[4];
	mp_idx_t *mi;

	fp = strcmp(fn, "-") == 0? stdin : fopen(fn, "rb");
	if (fp == 0) return 0;
	fread(magic, 1, 4, fp);
	if (strncmp(magic, MP_IDX_MAGIC, 4) != 0)
		return 0;
	mi = Kcalloc(0, mp_idx_t, 1);
	fread(&mi->opt, sizeof(mi->opt), 1, fp);
	ns_make_tables(mi->opt.trans_code);
	fread(&mi->n_kb, 8, 1, fp);
	mi->nt = mp_ntseq_restore(fp);
	mi->ki = Kmalloc(0, int64_t, mp_n_bucket(&mi->opt));
	mi->kb = Kmalloc(0, uint32_t, mi->n_kb);
	fread(mi->ki, 8, mp_n_bucket(&mi->opt), fp);
	fread(mi->kb, 4, mi->n_kb, fp);
	if (fp != stdin) fclose(fp);
	mi->bo = mp_idx_boff(mi->nt, mi->opt.bbit, &mi->n_block);
	if (mp_verbose >= 3)
		fprintf(stderr, "[M::%s@%.3f*%.2f] loaded the index\n", __func__, mp_realtime(), mp_percent_cpu());
	return mi;
}

mp_idx_t *mp_idx_load(const char *fn, const mp_idxopt_t *io, int32_t n_threads)
{
	int64_t is_idx;
	is_idx = mp_idx_is_idx(fn);
	if (is_idx < 0) return 0;
	return is_idx? mp_idx_restore(fn) : mp_idx_build(fn, io, n_threads);
}
