#ifndef MINIPROT_H
#define MINIPROT_H

#include <stdint.h>

#define MP_CODON_STD 0

#define MP_MAGIC "MPI\1"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct { uint64_t x, y; } mp128_t;
typedef struct { int32_t n, m; mp128_t *a; } mp128_v;
typedef struct { int32_t n, m; uint64_t *a; } mp64_v;

typedef struct {
	int32_t bbit; // block bit
	int32_t min_aa_len; // ignore ORFs shorter than this
	int32_t kmer, smer;
} mp_idxopt_t;

typedef struct {
	int64_t off, len;
	char *name;
} mp_ctg_t;

typedef struct {
	int32_t n_ctg, m_ctg;
	int32_t l_name;
	int64_t l_seq, m_seq;
	uint8_t *seq; // TODO: separate this into multiple blocks; low priority
	mp_ctg_t *ctg;
	char *name;
} mp_ntdb_t;

typedef struct {
	mp_idxopt_t opt;
	uint32_t n_block;
	mp_ntdb_t *nt;
	int64_t n_kb, *ki;
	uint32_t *bo, *kb;
} mp_idx_t;

extern int32_t mp_verbose;
extern char *mp_tab_nt_i2c, *mp_tab_aa_i2c;
extern uint8_t mp_tab_a2r[22], mp_tab_nt4[256], mp_tab_aa20[256], mp_tab_codon[64], mp_tab_codon13[64];

void mp_start(void);
void mp_make_tables(int codon_type);

void mp_idxopt_init(mp_idxopt_t *io);

int64_t mp_ntseq_get(const mp_ntdb_t *db, int32_t cid, int64_t st, int64_t en, int32_t rev, uint8_t *seq);
mp_idx_t *mp_idx_build(const char *fn, const mp_idxopt_t *io, int32_t n_threads);
mp_idx_t *mp_idx_read(const char *fn, const mp_idxopt_t *io, int32_t n_threads);
void mp_idx_destroy(mp_idx_t *mi);

int mp_idx_dump(const char *fn, const mp_idx_t *mi);
mp_idx_t *mp_idx_restore(const char *fn);

#ifdef __cplusplus
}
#endif

#endif
