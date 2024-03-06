#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "nasw.h"

char *ns_tab_nt_i2c = "ACGTN";
					// 01234

char *ns_tab_aa_i2c = "ARNDCQEGHILKMFPSTWYV*X";
					// 0123456789012345678901

uint8_t ns_tab_a2r[22] = { 0, 2, 4, 4, 6, 5, 5, 8, 3, 10, 11, 2, 11, 12, 7, 1, 1, 13, 12, 10, 14, 15 };
						// A  R  N  D  C  Q  E  G  H   I   L  K   M   F  P  S  T   W   Y   V   *   X

#define NS_MAX_TRANS_CODE 33
static const char *ns_tab_codon_all[NS_MAX_TRANS_CODE + 1] = {
   0,
 // 0123456789012345678901234567890123456789012345678901234567890123
 // A               C               G               T
 // A   C   G   T   A   C   G   T   A   C   G   T   A   C   G   T
 // ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLFX", // 1: The Standard Code; in order of AAA, AAC, AAG, AAT, ACA, ...
   "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLFX", // 2: The Vertebrate Mitochondrial Code
   "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLFX", // 3: The Yeast Mitochondrial Code
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLFX", // 4: The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
   "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLFX", // 5: The Invertebrate Mitochondrial Code
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLFX", // 6: Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear
   0, // 7
   0, // 8
   "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLFX", // 9: Echinoderm Mitochondrial; Flatworm Mitochondrial
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLFX", // 10: Euplotid Nuclear
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLFX", // 11: Bacterial, Archaeal and Plant Plastid
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLFX", // 12: Alternative Yeast Nuclear
   "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLFX", // 13: Ascidian Mitochondrial
   "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLFX", // 14: Alternative Flatworm Mitochondrial
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLFX", // 15: Blepharisma Macronuclear
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLFX", // 16: Chlorophycean Mitochondrial
   0, // 17
   0, // 18
   0, // 19
   0, // 20
   "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLFX", // 21: Trematode Mitochondrial
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLFX", // 22: Scenedesmus obliquus Mitochondrial
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLFX", // 23: Thraustochytrium Mitochondrial
   "KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLFX", // 24: Rhabdopleuridae Mitochondrial
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSGCWCLFLFX", // 25: Candidate Division SR1 and Gracilibacteria
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLALEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLFX", // 26: Pachysolen tannophilus Nuclear
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSSWCWCLFLFX", // 27: Karyorelict Nuclear
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSSWCWCLFLFX", // 28: Condylostoma Nuclear
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYYYYSSSS*CWCLFLFX", // 29: Mesodinium Nuclear
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVEYEYSSSS*CWCLFLFX", // 30: Peritrich Nuclear
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVEYEYSSSSWCWCLFLFX", // 31: Blastocrithidia Nuclear
   "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YWYSSSS*CWCLFLFX", // 32: Balanophoraceae Plastid
   "KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLFX"  // 33: Cephalodiscidae Mitochondrial
};

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

int ns_make_tables(int codon_type)
{
	const char *trans_tab;
	char *p;
	int i;
	if (codon_type < 0 || codon_type > NS_MAX_TRANS_CODE) return -1; // out of range
	trans_tab = ns_tab_codon_all[codon_type];
	if (trans_tab == 0) return -2; // not defined
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
		ns_tab_codon[i] = ns_tab_aa20[(uint8_t)trans_tab[i]];
		ns_tab_codon13[i] = ns_tab_a2r[ns_tab_codon[i]];
	}
	return 0;
}

/* See Sibley et al (2016)
 * Donor:
 *   P(GT.) = 53.58+45.10+0.37*0.73+0.06*0.33*0.5 = 98.96
 *   P(GTR) = 53.58+45.10*0.88+0.37+0.06*0.33*0.5*0.72 = 93.645128; s(GTY) = 2*log2((98.96-93.645128)/93.645128) = -8.278194172291354
 *   P(GC.) = 0.87+0.06*0.08*0.33 = 0.871584;  s(GC.) = 2*log2(0.871584/4/(93.645128/2)) = -15.494840831327044
 *   P(AT.) = 0.37*0.27+0.06*0.5*0.5 = 0.1149; s(AT.) = 2*log2(0.1149/4/(93.645128/2)) = -21.341362660546288
 *   P(...) = 0.053516; s(...) = 2*log2(0.053516/(64-4*3)/(93.645128/2)) = -30.9469153081545
 *   P(GGT.) = 53.58*0.7+45.10*0.92+0.37*0.73*0.18+0.06*0.33*0.5*0.55 = 79.05; s(HGT.) = 2*log2((98.96-79.05)/79.05) = -4.2840002505195445
 * Acceptor:
 *   P(.AG) = 64.55+29.01+5.98+0.37*0.73 = 99.8101
 *   P(YAG) = 64.55+29.01+0.37*0.73*0.9 = 93.80309; s(RAG) = 2*log2((99.8101-93.80309)/93.80309) = -7.929832954348969
 *   P(.AC) = 0.37*0.27 = 0.0999; s(.AC) = 2*log2(0.0999/4/(93.80309/2)) = -21.74987010896733
 *   P(...) = 0.09; s(...) = 2*log2(0.09/(64-4*2)/(93.80309/2)) = -29.6656993062333
 *
 */
void ns_opt_set_sp(ns_opt_t *opt, int32_t model)
{
	if (model == NS_S_MAMMAL) opt->sp[0] = 8, opt->sp[1] = 15, opt->sp[2] = 21, opt->sp[3] = 30, opt->sp[4] = 4, opt->sp[5] = 4;
	else if (model == NS_S_GENERIC) opt->sp[0] = 8, opt->sp[1] = 15, opt->sp[2] = 21, opt->sp[3] = 30, opt->sp[4] = opt->sp[5] = 0;
	else opt->sp[0] = opt->sp[1] = opt->sp[2] = opt->sp[3] = opt->sp[4] = opt->sp[5] = 0;
}

void ns_opt_init(ns_opt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	opt->go = 11, opt->ge = 1;
	opt->io = 29;
	opt->fs = 17;
	opt->xdrop = 100;
	opt->end_bonus = 5;
	ns_opt_set_sp(opt, NS_S_MAMMAL);
	opt->asize = 22;
	opt->ie_coef = .5f;
	opt->sc = ns_mat_blosum62;
	opt->nt4 = ns_tab_nt4;
	opt->aa20 = ns_tab_aa20;
	opt->codon = ns_tab_codon;
}

void ns_set_stop_sc(int32_t asize, int8_t *mat, int8_t pen)
{
	int32_t aa_stop = ns_tab_aa20['*'];
	int32_t i, score_ori = mat[aa_stop * asize + aa_stop];
	for (i = 0; i < asize; ++i)
		mat[aa_stop * asize + i] = mat[i * asize + aa_stop] = -pen;
	mat[aa_stop * asize + aa_stop] = score_ori;
}
