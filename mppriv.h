#ifndef MPPRIV_H
#define MPPRIV_H

#include <stdio.h>
#include "miniprot.h"
#include "nasw.h"
#include "kseq.h"

#define MP_DBG_NO_KALLOC   0x1
#define MP_DBG_QNAME       0x2
#define MP_DBG_NO_REFINE   0x4
#define MP_DBG_MORE_DP     0x8
#define MP_DBG_ANCHOR      0x10
#define MP_DBG_CHAIN       0x20

#define MP_PARENT_UNSET   (-1)
#define MP_PARENT_TMP_PRI (-2)

#ifndef kroundup64
#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern int8_t mp_mat_blosum62[484];

struct mp_bseq_file_s;
typedef struct mp_bseq_file_s mp_bseq_file_t;

typedef struct {
	int32_t l_seq;
	char *name, *seq, *comment;
} mp_bseq1_t;

// from bseq.c
mp_bseq_file_t *mp_bseq_open(const char *fn);
void mp_bseq_close(mp_bseq_file_t *fp);
mp_bseq1_t *mp_bseq_read(mp_bseq_file_t *fp, int64_t chunk_size, int with_comment, int *n_);

// from misc.c
char *mp_strdup(const char *src);
void radix_sort_mp64(uint64_t *st, uint64_t *en);
void radix_sort_mp128x(mp128_t *st, mp128_t *en);

// from sys.c
double mp_realtime(void);
double mp_cputime(void);
long mp_peakrss(void);
double mp_percent_cpu(void);

// from ntseq.c
mp_ntdb_t *mp_ntseq_read(const char *fn);
void mp_ntseq_destroy(mp_ntdb_t *db);
void mp_ntseq_dump(FILE *fp, const mp_ntdb_t *nt);
mp_ntdb_t *mp_ntseq_restore(FILE *fp);
int64_t mp_ntseq_get(const mp_ntdb_t *db, int32_t cid, int64_t st, int64_t en, int32_t rev, uint8_t *seq);
int64_t mp_ntseq_get_by_v(const mp_ntdb_t *nt, int32_t vid, int64_t st, int64_t en, uint8_t *seq);
int64_t mp_ntseq_spsc_get(const mp_ntdb_t *db, int32_t cid, int64_t st0, int64_t en0, int32_t rev, uint8_t *sc);
int64_t mp_ntseq_spsc_get_by_v(const mp_ntdb_t *nt, int32_t vid, int64_t st, int64_t en, uint8_t *ss);

// from sketch.c
void mp_sketch_nt4(void *km, const uint8_t *seq, int64_t len, int32_t min_aa_len, int32_t kmer, int32_t mod_bit, int32_t bbit, int64_t boff, mp64_v *a);
void mp_sketch_prot(void *km, const char *seq, int32_t len, int32_t kmer, int32_t mod_bit, mp64_v *a);

// from index.c
int32_t mp_idx_block2pos(const mp_idx_t *mi, uint32_t b);

// from hit.c
int32_t mp_cal_chn_sc_ungap(int32_t n_a, const uint64_t *a, int32_t kmer);
void mp_sort_reg(void *km, int *n_regs, mp_reg1_t *r);
void mp_set_parent(void *km, float mask_level, int mask_len, int n, mp_reg1_t *r, int sub_diff, int hard_mask_level);
void mp_select_sub(void *km, float pri_ratio, int min_diff, int best_n, int *n_, mp_reg1_t *r);
uint64_t *mp_cal_max_ext(void *km, int32_t n_reg, mp_reg1_t *reg, const uint64_t *a, int32_t min_ext, int32_t max_ext);
uint64_t *mp_collate_a(void *km, int32_t n_reg, mp_reg1_t *reg);
void mp_select_multi_exon(int32_t n, mp_reg1_t *r, int32_t single_penalty);

// from chain.c
uint64_t *mp_chain(int32_t max_dist_x, int32_t max_dist_y, int32_t bw, int32_t max_skip, int32_t max_iter, int32_t min_cnt, int32_t min_sc, float chn_coef_log,
				   int32_t is_spliced, int32_t kmer, int32_t bbit, int64_t n, uint64_t *a, int32_t *n_u_, uint64_t **_u, void *km);

mp_reg1_t *mp_reg_gen_from_block(void *km, const mp_idx_t *mi, int32_t n_u, const uint64_t *u, const uint64_t *a, int32_t *n_reg);

// from format.c
void mp_write_output(kstring_t *s, void *km, const mp_idx_t *mi, const mp_bseq1_t *seq, const mp_reg1_t *r, const mp_mapopt_t *opt, int64_t id, int32_t hit_idx);

// from align.c
void mp_align(void *km, const mp_mapopt_t *opt, const mp_idx_t *mi, int32_t len, const char *seq, mp_reg1_t *r, int32_t extl0, int32_t extr0);

static inline float mp_log2(float x) // NB: this doesn't work when x<2
{
	union { float f; uint32_t i; } z = { x };
	float log_2 = ((z.i >> 23) & 255) - 128;
	z.i &= ~(255 << 23);
	z.i += 127 << 23;
	log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
	return log_2;
}

static inline uint32_t mp_n_bucket(const mp_idxopt_t *io)
{
	return 1U << (io->kmer * MP_BITS_PER_AA - io->mod_bit);
}

#ifdef __cplusplus
}
#endif

#endif
