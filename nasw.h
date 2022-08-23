#ifndef NASW_H
#define NASW_H

#include <stdint.h>

#define NS_CIGAR_M	0
#define NS_CIGAR_I	1
#define NS_CIGAR_D	2
#define NS_CIGAR_N	3
#define NS_CIGAR_F	10 // 1bp frameshift
#define NS_CIGAR_G	11 // 2bp frameshift
#define NS_CIGAR_U	12
#define NS_CIGAR_V	13

#define NS_CIGAR_STR   "MIDNSHP=XBFGUV"

extern char *ns_tab_nt_i2c, *ns_tab_aa_i2c;
extern uint8_t ns_tab_a2r[22], ns_tab_nt4[256], ns_tab_aa20[256], ns_tab_aa13[256];
extern uint8_t ns_tab_codon[64], ns_tab_codon13[64];

typedef struct {
	int32_t go, ge, io, nc, fs;
	int32_t asize;
	int8_t *sc;
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

void ns_rst_init(ns_rst_t *r);
void ns_splice_s1(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, ns_rst_t *r);

#ifdef __cplusplus
}
#endif

#endif
