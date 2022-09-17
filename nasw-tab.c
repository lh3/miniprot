#include <string.h>
#include <ctype.h>
#include "nasw.h"

char *ns_tab_nt_i2c = "ACGTN";
					// 01234

char *ns_tab_aa_i2c = "ARNDCQEGHILKMFPSTWYV*X";
					// 0123456789012345678901

uint8_t ns_tab_a2r[22] = { 0, 2, 4, 4, 6, 5, 5, 8, 3, 10, 11, 2, 11, 12, 7, 1, 1, 13, 12, 10, 14, 15 };
						// A  R  N  D  C  Q  E  G  H   I   L  K   M   F  P  S  T   W   Y   V   *   X

char *ns_tab_codon_std = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLFX";
					   // 01234567890123456789012345678901234567890123456789012345678901234
					   // KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFFX <- this is the AGCT order

uint8_t ns_tab_nt4[256], ns_tab_aa20[256], ns_tab_aa13[256], ns_tab_codon[64], ns_tab_codon13[64];

int8_t ns_mat_blosum62[484] = { // 484 = 22*22
//	 A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  *  X
	 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-4, 0,
	-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-4,-1,
	-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,-4,-1,
	-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,-4,-1,
	 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-4,-2,
	-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,-4,-1,
	-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,-4,-1,
	 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-4,-1,
	-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,-4,-1,
	-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-4,-1,
	-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-1,
	-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,-4,-1,
	-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-4,-1,
	-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-4,-1,
	-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-4,-2,
	 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,-4, 0,
	 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-4, 0,
	-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-2,
	-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-4,-1,
	 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-4,-1,
	-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1,-4,
	 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-4,-1
};

void ns_make_tables(int codon_type)
{
	char *p;
	int i;
	memset(ns_tab_nt4, 4, 256);
	for (p = ns_tab_nt_i2c; *p; ++p)
		ns_tab_nt4[p - ns_tab_nt_i2c] = ns_tab_nt4[(uint8_t)toupper(*p)] = ns_tab_nt4[(uint8_t)tolower(*p)] = p - ns_tab_nt_i2c;
	memset(ns_tab_aa20, 21, 256);
	for (p = ns_tab_aa_i2c; *p; ++p)
		ns_tab_aa20[p - ns_tab_aa_i2c] = ns_tab_aa20[(uint8_t)toupper(*p)] = ns_tab_aa20[(uint8_t)tolower(*p)] = p - ns_tab_aa_i2c;
	memset(ns_tab_aa13, 15, 256);
	for (p = ns_tab_aa_i2c; *p; ++p)
		ns_tab_aa13[p - ns_tab_aa_i2c] = ns_tab_aa13[(uint8_t)toupper(*p)] = ns_tab_aa13[(uint8_t)tolower(*p)] = ns_tab_a2r[p - ns_tab_aa_i2c];
	for (i = 0; i < 64; ++i) {
		ns_tab_codon[i] = ns_tab_aa20[(uint8_t)ns_tab_codon_std[i]];
		ns_tab_codon13[i] = ns_tab_a2r[ns_tab_codon[i]];
	}
}

void ns_opt_init(ns_opt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	opt->go = 11, opt->ge = 1;
	opt->io = 31;
	opt->nc = 11;
	opt->fs = 17;
	opt->xdrop = 100;
	opt->end_bonus = 10;
	opt->asize = 22;
	opt->sc = ns_mat_blosum62;
	opt->nt4 = ns_tab_nt4;
	opt->aa20 = ns_tab_aa20;
	opt->codon = ns_tab_codon;
}
