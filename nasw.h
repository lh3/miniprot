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

typedef struct {
	int32_t n_cigar, m_cigar;
	int32_t score;
	uint32_t *cigar;
} ns_rst_t;

#endif
