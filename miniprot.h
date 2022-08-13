#ifndef MINIPROT_H
#define MINIPROT_H

#include <stdint.h>

#define MP_CODON_STD 0

#ifdef __cplusplus
extern "C" {
#endif

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

/*
 * table.c
 */
extern char *mp_tab_nt_i2c, *mp_tab_aa_i2c;
extern uint8_t mp_tab_a2r[22], mp_tab_nt4[256], mp_tab_aa20[256], mp_tab_codon[64];

void mp_make_tables(int codon_type);

/*
 * ntseq.c
 */
mp_ntdb_t *mp_ntseq_read(const char *fn);
int64_t mp_ntseq_get(const mp_ntdb_t *db, int32_t cid, int64_t st, int64_t en, int32_t rev, uint8_t *seq);
void mp_ntseq_get_orf(int64_t len, const uint8_t *seq);

#ifdef __cplusplus
}
#endif

#endif
