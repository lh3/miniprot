/*
  The MIT License

  Copyright (c) 2022-     Dana-Farber Cancer Institute

  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/
#ifndef NASW_H
#define NASW_H

#include <stdint.h>
#include <string.h>
#include "kalloc.h"

#define NS_CIGAR_M	0
#define NS_CIGAR_I	1
#define NS_CIGAR_D	2
#define NS_CIGAR_N	3  // phase-0 intron
// the following CIGAR operators are unique to nasw:
#define NS_CIGAR_F	10 // frameshift deletion
#define NS_CIGAR_G	11 // frameshift match
#define NS_CIGAR_U	12 // phase-1 intron
#define NS_CIGAR_V	13 // phase-2 intron
#define NS_CIGAR_E	14 // empty placeholder (not used)

#define NS_CIGAR_STR   "MIDNSHP=XBFGUVE"

#define NS_F_CIGAR      0x1
#define NS_F_EXT_LEFT   0x2
#define NS_F_EXT_RIGHT  0x4

#define NS_S_NONE       0
#define NS_S_GENERIC    1
#define NS_S_MAMMAL     2

extern char *ns_tab_nt_i2c, *ns_tab_aa_i2c;
extern uint8_t ns_tab_a2r[22], ns_tab_nt4[256], ns_tab_aa20[256], ns_tab_aa13[256];
extern uint8_t ns_tab_codon[64], ns_tab_codon13[64];
extern int8_t ns_mat_blosum62[484];

typedef struct {
	int32_t flag;
	int32_t go, ge, io, fs; // gap open, extension, intron open, frameshift
	int32_t xdrop, end_bonus; // xdrop for extension, and bonus for reaching end of proteins
	int32_t asize; // alphabet size; always 22 in the current implementation
	int32_t sp[6]; // 0:pos3, 1:GC-AC, 2:AT-AC, 3:other, 4:pos0, 5:poly-Y
	float ie_coef;
	const int8_t *sc; // 22x22 scoring matrix
	uint8_t *nt4, *aa20, *codon; // char-to-int convertion table and codon table
} ns_opt_t;

typedef struct {
	int32_t n_cigar, m_cigar; // number of cigar operators; max number
	int32_t nt_len, aa_len; // length of nt/aa sequence in the alignment
	int32_t score; // alignment score
	uint32_t *cigar; // CIGAR; allocated from the "km" memory block in ns_global_gs16()
} ns_rst_t;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Generate codon and conversion tables.
 *
 * This function should be called before invoking other nasw functions. It is
 * not thread safe.
 *
 * @param codon_type    NCBI genetic code
 *
 * @return 0 if tables were initialized; <0 if the translation table is not defined
 */
int ns_make_tables(int codon_type);

/**
 * Initialize alignment parameters
 *
 * @param opt      options to initialize
 */
void ns_opt_init(ns_opt_t *opt);

/**
 * Aligning a protein sequence to a nucleotide sequence with splicing and frameshift
 *
 * This function impelements striped dynamic programming (DP) using SSE2 on
 * x86_64 or NEON on ARM. It supports global alignment and left or right
 * extension. In the extension mode, only extension ends and alignment score
 * are computed. Users need to call the function again to get CIGAR.
 *
 * This function only works when the max score is not larger than 32767. If the
 * max score is larger, the 32-bit version ns_global_gs32() should be used
 * instead. ns_global_gs32() is twice as slow and uses twice as much memory and
 * it does not support extension alignment or end_bonus yet.
 *
 * The initial conditions of the DP are not technically correct. This issue may
 * affect gap openings at the beginning of the input sequences. The rest of
 * alignment should be ok.
 *
 * @param km        kalloc allocator; NULL to use the system malloc()
 * @param ns        the nucleotide sequence; can be normal char string or 0/1/2/3/4 encoding
 * @param nl        length of the nucleotide sequence
 * @param as        the protein sequence
 * @param al        the length of the protein sequence
 * @param opt       alignment settings and parameters
 * @param rst       (out) the result
 */
void ns_global_gs16(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, ns_rst_t *r);
void ns_global_gs32(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, ns_rst_t *r);

void ns_global_gs16b(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, const uint8_t *ss, ns_rst_t *r);
void ns_global_gs32b(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, const uint8_t *ss, ns_rst_t *r);

void ns_opt_set_sp(ns_opt_t *opt, int32_t model);

// Non-SIMD implementation with transposed DP matrix. For debugging only. Not intended for general developers.
void ns_splice_s1(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, ns_rst_t *r);

void ns_set_stop_sc(int32_t asize, int8_t *mat, int8_t score);

static inline uint32_t *ns_push_cigar(void *km, int32_t *n_cigar, int32_t *m_cigar, uint32_t *cigar, uint32_t op, int32_t len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf) || op == NS_CIGAR_F || op == NS_CIGAR_G) {
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
