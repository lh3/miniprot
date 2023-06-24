#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "kseq.h"
#include "mppriv.h"

static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

void mp_sprintf_lite(kstring_t *s, const char *fmt, ...) // FIXME: make it work for 64-bit integers
{
	char buf[16]; // for integer to string conversion
	const char *p, *q;
	va_list ap;
	va_start(ap, fmt);
	for (q = p = fmt; *p; ++p) {
		if (*p == '%') {
			if (p > q) str_copy(s, q, p);
			++p;
			if (*p == 'd') {
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 'u') {
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 's') {
				char *r = va_arg(ap, char*);
				str_copy(s, r, r + strlen(r));
			} else if (*p == 'c') {
				str_enlarge(s, 1);
				s->s[s->l++] = va_arg(ap, int);
			} else {
				fprintf(stderr, "ERROR: unrecognized type '%%%c'\n", *p);
				abort();
			}
			q = p + 1;
		}
	}
	if (p > q) str_copy(s, q, p);
	va_end(ap);
	s->s[s->l] = 0;
}

static void mp_write_cs(kstring_t *str, const mp_idx_t *mi, const char *aa, const mp_reg1_t *r)
{
	int32_t k, i, j, l, nl = 0, al = 0, l_tmp = 16;
	const mp_extra_t *e = r->p;
	int64_t l_nt;
	uint8_t *nt;
	char *tmp;
	if (e == 0) return;
	mp_sprintf_lite(str, "cs:Z:");
	for (k = 0; k < e->n_cigar; ++k) { // pre-calculate l_tmp
		int32_t op = e->cigar[k]&0xf, len = e->cigar[k]>>4, len3 = len * 3;
		if (op == NS_CIGAR_I || op == NS_CIGAR_F || op == NS_CIGAR_G) l_tmp = l_tmp > len? l_tmp : len;
		else if (op == NS_CIGAR_D) l_tmp = l_tmp > len3? l_tmp : len3;
	}
	tmp = Kmalloc(0, char, l_tmp + 16);
	nt = Kmalloc(0, uint8_t, r->ve - r->vs);
	l_nt = mp_ntseq_get_by_v(mi->nt, r->vid, r->vs, r->ve, nt);
	assert(l_nt == r->ve - r->vs);
	for (k = 0; k < e->n_cigar; ++k) {
		int32_t t, op = e->cigar[k]&0xf, len = e->cigar[k]>>4, len3 = len * 3;
		if (op == NS_CIGAR_M) {
			for (i = nl, j = al, l = 0, t = 0; l < len; ++l, ++j, i += 3) {
				uint8_t nt_aa, aa_aa, codon = nt[i]<<4 | nt[i+1]<<2 | nt[i+2];
				nt_aa = nt[i] > 3 || nt[i+1] > 3 || nt[i+2] > 3? ns_tab_aa20['X'] : ns_tab_codon[codon];
				aa_aa = ns_tab_aa20[(uint8_t)aa[j]];
				if (nt_aa != aa_aa) {
					if (t > 0) mp_sprintf_lite(str, ":%d", t);
					tmp[0] = "acgtn"[nt[i]], tmp[1] = "acgtn"[nt[i+1]], tmp[2] = "acgtn"[nt[i+2]], tmp[3] = toupper(aa[j]), tmp[4] = 0;
					mp_sprintf_lite(str, "*%s", tmp);
					t = 0;
				} else ++t;
			}
			if (t > 0) mp_sprintf_lite(str, ":%d", t); // t is the length of exact match
			nl += len3, al += len;
		} else if (op == NS_CIGAR_I) {
			for (j = t = 0; j < len; ++j)
				tmp[t++] = toupper(aa[al + j]);
			tmp[t] = 0;
			mp_sprintf_lite(str, "+%s", tmp);
			al += len;
		} else if (op == NS_CIGAR_D) {
			for (i = t = 0; i < len3; ++i)
				tmp[t++] = "acgtn"[nt[nl + i]];
			tmp[t] = 0;
			mp_sprintf_lite(str, "-%s", tmp);
			nl += len3;
		} else if (op == NS_CIGAR_F) {
			for (i = t = 0; i < len; ++i)
				tmp[t++] = "acgtn"[nt[nl + i]];
			tmp[t] = 0;
			mp_sprintf_lite(str, "-%s", tmp);
			nl += len;
		} else if (op == NS_CIGAR_G) {
			for (i = t = 0; i < len; ++i)
				tmp[t++] = "acgtn"[nt[nl + i]];
			tmp[t++] = toupper(aa[al]);
			tmp[t] = 0;
			mp_sprintf_lite(str, "*%s", tmp);
			nl += len, ++al;
		} else if (op == NS_CIGAR_N || op == NS_CIGAR_U || op == NS_CIGAR_V) {
			int32_t lshift = op == NS_CIGAR_N? 0 : op == NS_CIGAR_U? 1 : 2;
			int32_t rshift = lshift == 0? 0 : 3 - lshift;
			if (lshift > 0) {
				for (i = t = 0; i < lshift; ++i)
					tmp[t++] = "acgtn"[nt[nl + i]];
				tmp[t++] = toupper(aa[al]);
				tmp[t] = 0;
				mp_sprintf_lite(str, "*%s", tmp);
			}
			mp_sprintf_lite(str, "~%c%c%d%c%c", "acgtn"[nt[nl + lshift]], "acgtn"[nt[nl + lshift + 1]], len - (lshift + rshift),
							"acgtn"[nt[nl + len - rshift - 2]], "acgtn"[nt[nl + len - rshift - 1]]);
			if (rshift > 0) {
				for (i = t = 0; i < rshift; ++i)
					tmp[t++] = "acgtn"[nt[nl + len - rshift + i]];
				tmp[t] = 0;
				mp_sprintf_lite(str, "-%s", tmp);
			}
			if (lshift) ++al;
			nl += len;
		}
	}
	assert(nl == r->ve - r->vs);
	assert(al == r->qe - r->qs);
	kfree(0, nt);
	kfree(0, tmp);
}

static void mp_write_residue(kstring_t *out, const mp_idx_t *mi, const mp_mapopt_t *opt, const char *aa, const mp_reg1_t *r)
{
	const mp_extra_t *e = r->p;
	int32_t k, w, w4, al = r->qs, nl = 0, l_tmp = 8; // enough to keep "##ATN\t\0"
	int32_t max_flank = opt->max_intron_flank;
	int64_t l_nt;
	uint8_t *nt;
	char *str[5];
	if (e == 0) return;
	for (k = 0; k < e->n_cigar; ++k) { // estimate the string length
		int32_t op = e->cigar[k]&0xf, len = e->cigar[k]>>4, len3 = len * 3;
		if (op == NS_CIGAR_M || op == NS_CIGAR_I || op == NS_CIGAR_D) l_tmp += len3;
		else if (op == NS_CIGAR_F || op == NS_CIGAR_G) l_tmp += len;
		else if (op == NS_CIGAR_N || op == NS_CIGAR_U || op == NS_CIGAR_V) {
			int32_t intron_len = op == NS_CIGAR_N? len : len - 3;
			l_tmp += intron_len <= max_flank * 2? len : max_flank * 2 + 20;
		}
	}
	l_tmp += 3;
	str[0] = Kmalloc(0, char, 5 * l_tmp);
	str[1] = str[0] + l_tmp;
	str[2] = str[1] + l_tmp;
	str[3] = str[2] + l_tmp;
	str[4] = str[3] + l_tmp;
	strcpy(str[0], "##ATN\t");
	strcpy(str[1], "##ATA\t");
	strcpy(str[2], "##AAS\t");
	strcpy(str[3], "##AQA\t");
	strcpy(str[4], "##STA\t");
	nt = Kmalloc(0, uint8_t, r->ve - r->vs + 3);
	l_nt = mp_ntseq_get_by_v(mi->nt, r->vid, r->vs, r->ve + 3, nt);
	assert(l_nt >= r->ve - r->vs);
	for (k = 0, w = w4 = 6; k < e->n_cigar; ++k) {
		int32_t i, j, l, op = e->cigar[k]&0xf, len = e->cigar[k]>>4, len3 = len * 3;
		assert(w < l_tmp);
		if (op == NS_CIGAR_M) {
			for (i = nl, j = al, l = 0; l < len; ++l, ++j, i += 3, w += 3) {
				uint8_t nt_aa, aa_aa, codon = nt[i]<<4 | nt[i+1]<<2 | nt[i+2];
				int32_t s;
				nt_aa = nt[i] > 3 || nt[i+1] > 3 || nt[i+2] > 3? ns_tab_aa20['X'] : ns_tab_codon[codon];
				aa_aa = ns_tab_aa20[(uint8_t)aa[j]];
				s = opt->mat[nt_aa * opt->asize + aa_aa];
				str[0][w] = "ACGTN"[nt[i]], str[0][w+1] = "ACGTN"[nt[i+1]], str[0][w+2] = "ACGTN"[nt[i+2]];
				str[1][w] = str[4][w4++] = ns_tab_aa_i2c[nt_aa], str[1][w+1] = str[1][w+2] = '.';
				str[2][w] = nt_aa == aa_aa? '|' : s > 0? '+' : ' ', str[2][w+1] = ' ', str[2][w+2] = ' ';
				str[3][w] = toupper(aa[j]), str[3][w+1] = str[3][w+2] = ' ';
			}
			nl += len3, al += len;
		} else if (op == NS_CIGAR_I) {
			for (j = 0; j < len; ++j, w += 3) {
				str[0][w] = str[0][w+1] = str[0][w+2] = '-';
				str[1][w] = '-', str[1][w+1] = str[1][w+2] = '.';
				str[2][w] = str[2][w+1] = str[2][w+2] = ' ';
				str[3][w] = toupper(aa[al + j]), str[3][w+1] = str[3][w+2] = ' ';
			}
			al += len;
		} else if (op == NS_CIGAR_D) {
			for (l = 0, i = nl; l < len; ++l, i += 3, w += 3) {
				uint8_t nt_aa, codon = nt[i]<<4 | nt[i+1]<<2 | nt[i+2];
				nt_aa = nt[i] > 3 || nt[i+1] > 3 || nt[i+2] > 3? ns_tab_aa20['X'] : ns_tab_codon[codon];
				str[0][w] = "ACGTN"[nt[i]], str[0][w+1] = "ACGTN"[nt[i+1]], str[0][w+2] = "ACGTN"[nt[i+2]];
				str[1][w] = str[4][w4++] = ns_tab_aa_i2c[nt_aa], str[1][w+1] = str[1][w+2] = '.';
				str[2][w] = str[2][w+1] = str[2][w+2] = ' ';
				str[3][w] = '-', str[3][w+1] = str[3][w+2] = ' ';
			}
			nl += len3;
		} else if (op == NS_CIGAR_F) {
			for (l = 0, i = nl; l < len; ++l, ++i, ++w)
				str[0][w] = "ACGTN"[nt[i]], str[1][w] = '!', str[2][w] = str[3][w] = ' ';
			nl += len;
		} else if (op == NS_CIGAR_G) {
			for (l = 0, i = nl; l < len; ++l, ++i, ++w) {
				str[0][w] = "ACGTN"[nt[i]], str[1][w] = '$', str[2][w] = ' ';
				str[3][w] = l == 0? toupper(aa[al]) : ' ';
			}
			nl += len, ++al;
		} else if (op == NS_CIGAR_N || op == NS_CIGAR_U || op == NS_CIGAR_V) {
			int32_t intron_len = op == NS_CIGAR_N? len : len - 3;
			if (op == NS_CIGAR_U || op == NS_CIGAR_V) { // phase-1 or phase-2 intron; print the first half of split codon
				uint8_t n1, n2, n3, codon, nt_aa, aa_aa;
				int32_t s;
				if (op == NS_CIGAR_U) n1 = nt[nl], n2 = nt[nl + len - 2], n3 = nt[nl + len - 1];
				else n1 = nt[nl], n2 = nt[nl + 1], n3 = nt[nl + len - 1];
				codon = n1<<4 | n2<<2 | n3;
				nt_aa = n1 > 3 || n2 > 3 || n3 > 3? ns_tab_aa20['X'] : ns_tab_codon[codon];
				aa_aa = ns_tab_aa20[(uint8_t)aa[al]];
				s = opt->mat[nt_aa * opt->asize + aa_aa];
				str[0][w] = "ACGTN"[nt[nl]];
				str[1][w] = str[4][w4++] = ns_tab_aa_i2c[nt_aa];
				str[2][w] = nt_aa == aa_aa? '|' : s > 0? '+' : ' ';
				str[3][w] = toupper(aa[al]);
				++w, ++nl;
				if (op == NS_CIGAR_V) {
					str[0][w] = "ACGTN"[nt[nl]], str[1][w] = '.', str[2][w] = str[3][w] = ' ';
					++w, ++nl;
				}
				++al;
			}
			if (intron_len <= max_flank * 2) {
				for (l = 0, i = nl; l < intron_len; ++l, ++i, ++w)
					str[0][w] = "acgtn"[nt[i]], str[1][w] = str[2][w] = str[3][w] = ' ';
			} else {
				int32_t il;
				for (l = 0, i = nl; l < max_flank; ++l, ++i, ++w)
					str[0][w] = "acgtn"[nt[i]], str[1][w] = str[2][w] = str[3][w] = ' ';
				str[0][w] = '~', str[1][w] = str[2][w] = ' ', str[3][w++] = ' ';
				il = snprintf(&str[0][w], l_tmp - w, "%d", intron_len);
				for (l = 0; l < il; ++l, ++w)
					str[1][w] = str[2][w] = str[3][w] = ' ';
				str[0][w] = '~', str[1][w] = str[2][w] = ' ', str[3][w++] = ' ';
				for (l = 0, i = nl + intron_len - max_flank; l < max_flank; ++l, ++i, ++w)
					str[0][w] = "acgtn"[nt[i]], str[1][w] = str[2][w] = str[3][w] = ' ';
			}
			nl += intron_len;
			if (op == NS_CIGAR_U || op == NS_CIGAR_V) { // phase-1 or phase-2 intron; print the second half of split codon
				str[0][w] = "ACGTN"[nt[nl]], str[1][w] = '.', str[2][w] = str[3][w] = ' ';
				++w, ++nl;
				if (op == NS_CIGAR_U) {
					str[0][w] = "ACGTN"[nt[nl]], str[1][w] = '.', str[2][w] = str[3][w] = ' ';
					++w, ++nl;
				}
			}
		}
	}
	assert(nl == r->ve - r->vs);
	if (l_nt == r->ve - r->vs + 3) { // then print one more codon
		uint8_t n1 = nt[nl], n2 = nt[nl + 1], n3 = nt[nl + 2], nt_aa, codon;
		codon = n1<<4 | n2<<2 | n3;
		nt_aa = n1 > 3 || n2 > 3 || n3 > 3? ns_tab_aa20['X'] : ns_tab_codon[codon];
		str[0][w] = "ACGTN"[n1], str[0][w+1] = "ACGTN"[n2], str[0][w+2] = "ACGTN"[n3];
		str[1][w] = str[4][w4++] = ns_tab_aa_i2c[nt_aa], str[1][w+1] = str[1][w+2] = '.';
		str[2][w] = str[2][w+1] = str[2][w+2] = ' ';
		str[3][w] = str[3][w+1] = str[3][w+2] = ' ';
		w += 3;
	}
	str[0][w] = str[1][w] = str[2][w] = str[3][w] = 0, str[4][w4] = 0;
	kfree(0, nt);
	if (opt->flag & MP_F_SHOW_RESIDUE)
		mp_sprintf_lite(out, "%s\n%s\n%s\n%s\n", str[0], str[1], str[2], str[3]);
	if (opt->flag & MP_F_SHOW_TRANS)
		mp_sprintf_lite(out, "%s\n", str[4]);
	kfree(0, str[0]);
}

void mp_write_paf(kstring_t *s, const mp_idx_t *mi, const mp_bseq1_t *seq, const mp_reg1_t *r, int32_t gff_out)
{
	const mp_ctg_t *ctg;
	if (gff_out) mp_sprintf_lite(s, "##PAF\t");
	if (r == 0) {
		mp_sprintf_lite(s, "%s\t%d\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\n", seq->name, seq->l_seq);
		return;
	}
	ctg = &mi->nt->ctg[r->vid>>1];
	mp_sprintf_lite(s, "%s\t%d\t%d\t%d\t%c\t%s\t%d\t", seq->name, seq->l_seq, r->qs, r->qe, "+-"[r->vid&1], ctg->name, (int)ctg->len);
	if (r->vid&1) mp_sprintf_lite(s, "%d\t%d\t", (int)(ctg->len - r->ve), (int)(ctg->len - r->vs)); // FIXME: make it work for 64-bit integers
	else mp_sprintf_lite(s, "%d\t%d\t", (int)r->vs, (int)r->ve);
	if (r->p) {
		int32_t k;
		mp_sprintf_lite(s, "%d\t%d\t0\tAS:i:%d\tms:i:%d\tnp:i:%d\tfs:i:%d\tst:i:%d\tda:i:%d\tdo:i:%d\t",
			r->p->n_iden * 3, r->p->blen, r->p->dp_score, r->p->dp_max, r->p->n_plus, r->p->n_fs, r->p->n_stop, r->p->dist_start, r->p->dist_stop);
		mp_sprintf_lite(s, "cg:Z:");
		for (k = 0; k < r->p->n_cigar; ++k)
			mp_sprintf_lite(s, "%d%c", r->p->cigar[k]>>4, NS_CIGAR_STR[r->p->cigar[k]&0xf]);
	} else mp_sprintf_lite(s, "%d\t%d\t%d", r->chn_sc, r->chn_sc_ungap, r->cnt);
	mp_sprintf_lite(s, "\t");
	mp_write_cs(s, mi, &seq->seq[r->qs], r);
	mp_sprintf_lite(s, "\n");
}

static void mp_write_gff(kstring_t *s, void *km, const mp_idx_t *mi, const mp_bseq1_t *seq, const mp_reg1_t *r, const char *gff_prefix, int64_t id, const char *qname, int32_t gff_delim, int32_t hit_idx)
{
	const mp_ctg_t *ctg;
	const mp_feat_t *f, *feat = r->feat;
	int64_t vs, ve, ve_mRNA;
	int32_t j, has_stop = 0;
	char dec[16], *id_str = 0;

	if (r == 0 || r->p == 0) return;
	has_stop = (r->qe == seq->l_seq && r->p->dist_stop == 0);
	ve_mRNA = has_stop? r->ve + 3 : r->ve;

	if (gff_delim >= 33 && gff_delim <= 126 && hit_idx >= 0) {
		id_str = Kmalloc(km, char, strlen(qname) + 16);
		sprintf(id_str, "%s%c%d", qname, gff_delim, hit_idx);
	} else {
		id_str = Kmalloc(km, char, strlen(gff_prefix) + 16);
		sprintf(id_str, "%s%.6ld", gff_prefix, (long)id);
	}
	ctg = &mi->nt->ctg[r->vid>>1];
	vs = r->vid&1? ctg->len - ve_mRNA : r->vs;
	ve = r->vid&1? ctg->len - r->vs   : ve_mRNA;
	mp_sprintf_lite(s, "%s\tminiprot\tmRNA\t%d\t%d\t%d\t%c\t.\tID=%s;Rank=%d", ctg->name, (int)vs + 1, (int)ve, r->p->dp_max, "+-"[r->vid&1], id_str, hit_idx);
	snprintf(dec, 16, "%.4f", (double)r->p->n_iden * 3 / r->p->blen);
	mp_sprintf_lite(s, ";Identity=%s", dec);
	snprintf(dec, 16, "%.4f", (double)r->p->n_plus * 3 / r->p->blen);
	mp_sprintf_lite(s, ";Positive=%s", dec);
	if (r->p->n_fs > 0) mp_sprintf_lite(s, ";Frameshift=%d", r->p->n_fs);
	if (r->p->n_stop > 0) mp_sprintf_lite(s, ";StopCodon=%d", r->p->n_stop);
	mp_sprintf_lite(s, ";Target=%s %d %d\n", seq->name, r->qs + 1, r->qe);

	for (j = 0; j < r->n_feat; ++j) {
		f = &feat[j];
		vs = r->vid&1? ctg->len - f->ve : f->vs;
		ve = r->vid&1? ctg->len - f->vs : f->ve;
		mp_sprintf_lite(s, "%s\tminiprot\t%s\t%d\t%d\t%d\t%c\t%d\tParent=%s;Rank=%d", ctg->name, f->type == MP_FEAT_STOP? "stop_codon" : "CDS",
			(int)vs + 1, (int)ve, f->score, "+-"[r->vid&1], f->phase, id_str, hit_idx);
		if (f->type == MP_FEAT_CDS) {
			snprintf(dec, 16, "%.4f", (double)f->n_iden * 3 / f->blen);
			mp_sprintf_lite(s, ";Identity=%s", dec);
			if (f->acceptor[0] && strncmp(f->acceptor, "AG", 2) != 0) mp_sprintf_lite(s, ";Acceptor=%c%c", f->acceptor[0], f->acceptor[1]);
			if (f->donor[0]    && strncmp(f->donor,    "GT", 2) != 0) mp_sprintf_lite(s, ";Donor=%c%c",    f->donor[0],    f->donor[1]);
			if (f->n_fs > 0) mp_sprintf_lite(s, ";Frameshift=%d", f->n_fs);
			if (f->n_stop > 0) mp_sprintf_lite(s, ";StopCodon=%d", f->n_stop);
			mp_sprintf_lite(s, ";Target=%s %d %d", seq->name, f->qs + 1, f->qe);
		}
		mp_sprintf_lite(s, "\n");
	}
	kfree(km, id_str);
}

static void mp_write_gtf(kstring_t *s, void *km, const mp_idx_t *mi, const mp_bseq1_t *seq, const mp_reg1_t *r, const char *gff_prefix, int64_t id, const char *qname)
{
	const mp_ctg_t *ctg;
	const mp_feat_t *f, *feat = r->feat;
	int64_t vs, ve, ve_mRNA;
	int32_t j, has_stop = 0;
	char *id_g_str = 0, *id_t_str = 0;

	if (r == 0 || r->p == 0) return;
	has_stop = (r->qe == seq->l_seq && r->p->dist_stop == 0);
	ve_mRNA = has_stop? r->ve + 3 : r->ve;

	id_g_str = Kmalloc(km, char, strlen(gff_prefix) + 16);
	id_t_str = Kmalloc(km, char, strlen(gff_prefix) + 16);
	sprintf(id_g_str, "%sG%.6ld", gff_prefix, (long)id);
	sprintf(id_t_str, "%sT%.6ld", gff_prefix, (long)id);
	ctg = &mi->nt->ctg[r->vid>>1];
	vs = r->vid&1? ctg->len - ve_mRNA : r->vs;
	ve = r->vid&1? ctg->len - r->vs   : ve_mRNA;
	mp_sprintf_lite(s, "%s\tminiprot\tgene\t%d\t%d\t%d\t%c\t.\tgene_id \"%s\";\n", ctg->name, (int)vs + 1, (int)ve, r->p->dp_max, "+-"[r->vid&1], id_g_str);
	mp_sprintf_lite(s, "%s\tminiprot\ttranscript\t%d\t%d\t%d\t%c\t.\ttranscript_id \"%s\"; gene_id \"%s\";\n", ctg->name, (int)vs + 1, (int)ve, r->p->dp_max, "+-"[r->vid&1], id_t_str, id_g_str);

	for (j = 0; j < r->n_feat; ++j) {
		int64_t vs2, ve2;
		f = &feat[j];
		if (f->type != MP_FEAT_CDS) continue;
		vs2 = vs = r->vid&1? ctg->len - f->ve : f->vs;
		ve2 = ve = r->vid&1? ctg->len - f->vs : f->ve;
		if (f->ve == r->ve) { // last exon; then adjust for stop codon
			if (r->vid&1) vs2 = ctg->len - ve_mRNA;
			else ve2 = ve_mRNA;
		}
		mp_sprintf_lite(s, "%s\tminiprot\texon\t%d\t%d\t%d\t%c\t.\ttranscript_id \"%s\"; gene_id \"%s\";\n", ctg->name, (int)vs2 + 1, (int)ve2, f->score, "+-"[r->vid&1], id_t_str, id_g_str);
		mp_sprintf_lite(s, "%s\tminiprot\tCDS\t%d\t%d\t%d\t%c\t%d\ttranscript_id \"%s\"; gene_id \"%s\";\n", ctg->name, (int)vs + 1,  (int)ve,  f->score, "+-"[r->vid&1], f->phase, id_t_str, id_g_str);
	}
	kfree(km, id_g_str);
	kfree(km, id_t_str);
}

void mp_write_output(kstring_t *s, void *km, const mp_idx_t *mi, const mp_bseq1_t *seq, const mp_reg1_t *r, const mp_mapopt_t *opt, int64_t id, int32_t hit_idx)
{
	s->l = 0;
	if (r == 0) {
		if (opt->flag&MP_F_SHOW_UNMAP)
			mp_write_paf(s, mi, seq, 0, opt->flag&MP_F_GFF);
	} else if (opt->flag&MP_F_GTF) {
		if (opt->flag & (MP_F_SHOW_RESIDUE|MP_F_SHOW_TRANS)) {
			mp_write_paf(s, mi, seq, r, opt->flag&MP_F_GTF);
			mp_write_residue(s, mi, opt, seq->seq, r);
		}
		mp_write_gtf(s, km, mi, seq, r, opt->gff_prefix, id, seq->name);
	} else {
		if (!(opt->flag&MP_F_NO_PAF))
			mp_write_paf(s, mi, seq, r, opt->flag&MP_F_GFF);
		if (opt->flag & (MP_F_SHOW_RESIDUE|MP_F_SHOW_TRANS))
			mp_write_residue(s, mi, opt, seq->seq, r);
		if (opt->flag&MP_F_GFF)
			mp_write_gff(s, km, mi, seq, r, opt->gff_prefix, id, seq->name, opt->gff_delim, hit_idx);
	}
}
