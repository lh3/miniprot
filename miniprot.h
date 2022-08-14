#ifndef MINIPROT_H
#define MINIPROT_H

#include <stdint.h>

#define MP_CODON_STD 0

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
	int64_t l_seq, m_seq;
	uint8_t *seq; // TODO: separate this into multiple blocks; low priority
	mp_ctg_t *ctg;
} mp_ntdb_t;

typedef struct {
	mp_ntdb_t *nt;
	uint32_t *boff;
} mp_idx_t;

/*
 * options.c
 */
void mp_idxopt_init(mp_idxopt_t *io);

/*
 * table.c
 */
extern char *mp_tab_nt_i2c, *mp_tab_aa_i2c;
extern uint8_t mp_tab_a2r[22], mp_tab_nt4[256], mp_tab_aa20[256], mp_tab_codon[64], mp_tab_codon13[64];

void mp_make_tables(int codon_type);

/*
 * ntseq.c
 */
int64_t mp_ntseq_get(const mp_ntdb_t *db, int32_t cid, int64_t st, int64_t en, int32_t rev, uint8_t *seq);
void mp_idx_proc_seq(void *km, int64_t len, const uint8_t *seq, int32_t min_aa_len, int32_t kmer, int32_t smer, int32_t bbit, int64_t boff, mp64_v *a);
mp_idx_t *mp_index(const mp_idxopt_t *io, const char *fn);

#ifdef __cplusplus
}
#endif

#endif
