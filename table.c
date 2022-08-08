#include <string.h>
#include <ctype.h>
#include "miniprot.h"

char *mp_tab_nt_i2c = "ACGTN";
					// 01234

char *mp_tab_aa_i2c = "ARNDCQEGHILKMFPSTWYV*X";
                    // 0123456789012345678901

uint8_t mp_tab_a2r[22] = { 0, 2, 4, 4, 6, 5, 5, 8, 3, 10, 11, 2, 11, 12, 7, 1, 1, 13, 12, 10, 14, 15 };
                      // A  R  N  D  C  Q  E  G  H   I   L  K   M   F  P  S  T   W   Y   V   *   X

char *mp_tab_codon_std = "KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFFX";
					   // 01234567890123456789012345678901234567890123456789012345678901234

uint8_t mp_tab_nt4[256], mp_tab_aa20[256];

void mp_make_tables(int codon_type)
{
	char *p;
	memset(mp_tab_nt4, 4, 256);
	for (p = mp_tab_nt_i2c; *p; ++p)
		mp_tab_nt4[(uint8_t)toupper(*p)] = mp_tab_nt4[(uint8_t)tolower(*p)] = p - mp_tab_nt_i2c;
	memset(mp_tab_aa20, 21, 256);
	for (p = mp_tab_aa_i2c; *p; ++p)
		mp_tab_aa20[(uint8_t)toupper(*p)] = mp_tab_aa20[(uint8_t)tolower(*p)] = p - mp_tab_aa_i2c;
}
