#ifndef NASW_H
#define NASW_H

#include <stdint.h>
#include <string.h>
#include "kalloc.h"

#define NS_CIGAR_M	0
#define NS_CIGAR_I	1
#define NS_CIGAR_D	2
#define NS_CIGAR_N	3  // phase-0 intron
#define NS_CIGAR_F	10 // frameshift
#define NS_CIGAR_G	11 // frameshift
#define NS_CIGAR_U	12 // phase-1 intron (TODO: may change)
#define NS_CIGAR_V	13 // phase-2 intron

#define NS_CIGAR_STR   "MIDNSHP=XBFGUV"

#define NS_F_CIGAR    0x1

extern char *ns_tab_nt_i2c, *ns_tab_aa_i2c;
extern uint8_t ns_tab_a2r[22], ns_tab_nt4[256], ns_tab_aa20[256], ns_tab_aa13[256];
extern uint8_t ns_tab_codon[64], ns_tab_codon13[64];
extern int8_t ns_mat_blosum62[484];

typedef struct {
	int32_t flag;
	int32_t go, ge, io, nc, fs;
	int32_t asize;
	const int8_t *sc;
	uint8_t *nt4, *aa20, *codon;
} ns_opt_t;

typedef struct {
	int32_t n_cigar, m_cigar;
	int32_t score;
	uint32_t *cigar;
} ns_rst_t;

#ifdef __cplusplus
extern "C" {
#endif

void ns_make_tables(int codon_type);
void ns_opt_init(ns_opt_t *opt);

void ns_splice_s1(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, ns_rst_t *r);
void ns_splice_i16(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, ns_rst_t *r);

void ns_global_gs16(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, ns_rst_t *r);
void ns_global_gs32(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, ns_rst_t *r);

static inline uint32_t *ns_push_cigar(void *km, int32_t *n_cigar, int32_t *m_cigar, uint32_t *cigar, uint32_t op, int32_t len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			(*m_cigar) += ((*m_cigar)>>1) + 8;
			cigar = Krealloc(km, uint32_t, cigar, *m_cigar);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

static inline void ns_rst_init(ns_rst_t *r)
{
	memset(r, 0, sizeof(*r));
}

#ifdef __cplusplus
}
#endif

#endif
