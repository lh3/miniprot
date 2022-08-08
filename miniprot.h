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
extern uint8_t mp_tab_a2r[22], mp_tab_nt4[256], mp_tab_aa20[256];

void mp_make_tables(int codon_type);

/*
 * ntseq.c
 */
mp_ntdb_t *mp_ntseq_read(const char *fn);

#ifdef __cplusplus
}
#endif

#endif
