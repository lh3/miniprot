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
			if (t > 0) mp_sprintf_lite(str, ":%d", t);
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
		mp_sprintf_lite(s, "%d\t%d\t0\tAS:i:%d\tas:i:%d\tnp:i:%d\tda:i:%d\tdo:i:%d\t", r->p->n_iden * 3, r->p->clen, r->p->dp_max, r->p->aa_score, r->p->n_plus, r->p->dist_start, r->p->dist_stop);
		mp_sprintf_lite(s, "cg:Z:");
		for (k = 0; k < r->p->n_cigar; ++k)
			mp_sprintf_lite(s, "%d%c", r->p->cigar[k]>>4, NS_CIGAR_STR[r->p->cigar[k]&0xf]);
	} else mp_sprintf_lite(s, "%d\t%d", r->chn_sc, r->cnt);
	mp_sprintf_lite(s, "\t");
	mp_write_cs(s, mi, &seq->seq[r->qs], r);
	mp_sprintf_lite(s, "\n");
}

#define MP_FEAT_EXON  0
#define MP_FEAT_STOP  1

typedef struct {
	int64_t vs, ve;
	int32_t qs, qe;
	int32_t phase;
	int32_t n_fs;
	int32_t type;
} mp_feat_t;

void mp_write_gff(kstring_t *s, void *km, const mp_idx_t *mi, const mp_bseq1_t *seq, const mp_reg1_t *r, const char *gff_prefix, int64_t id, const char *qname)
{
	const mp_ctg_t *ctg;
	mp_feat_t *feat, *f;
	int32_t k, j, n_intron = 0, n_feat;
	int64_t vs, ve, ve_mRNA;
	int32_t qs, qe, phase, n_fs, tot_fs = 0, has_stop = 0;
	char buf[16], dec[16];

	if (r == 0 || r->p == 0) return;
	has_stop = (r->qe == seq->l_seq && r->p->dist_stop == 0);
	ve_mRNA = has_stop? r->ve + 3 : r->ve;
	for (k = 0; k < r->p->n_cigar; ++k) {
		int32_t op = r->p->cigar[k]&0xf;
		if (op == NS_CIGAR_N || op == NS_CIGAR_U || op == NS_CIGAR_V)
			++n_intron;
	}
	n_feat = n_intron + 1 + (!!has_stop);
	feat = Kmalloc(km, mp_feat_t, n_feat);

	vs = ve = r->vs, qs = qe = r->qs, phase = 0, n_fs = 0;
	for (k = j = 0; k < r->p->n_cigar; ++k) {
		int32_t op = r->p->cigar[k]&0xf, len = r->p->cigar[k]>>4, len3 = len * 3;
		if (op == NS_CIGAR_M) {
			ve += len3, qe += len;
		} else if (op == NS_CIGAR_I) {
			qe += len;
		} else if (op == NS_CIGAR_D) {
			ve += len3;
		} else if (op == NS_CIGAR_F) {
			ve += len, n_fs += len, tot_fs += len;
		} else if (op == NS_CIGAR_G) {
			ve += len, ++qe, n_fs += len, tot_fs += len;
		} else if (op == NS_CIGAR_N || op == NS_CIGAR_U || op == NS_CIGAR_V) {
			f = &feat[j++];
			f->phase = phase;
			f->vs = vs, f->qs = qs, f->qe = qe, f->n_fs = n_fs, f->type = MP_FEAT_EXON;
			if (op == NS_CIGAR_N) {
				f->ve = ve;
				vs = ve + len, phase = 0;
			} else if (op == NS_CIGAR_U) {
				f->ve = ve + 1;
				vs = ve + len - 2, phase = 2;
			} else if (op == NS_CIGAR_V) {
				f->ve = ve + 2;
				vs = ve + len - 1, phase = 1;
			}
			qs = qe, n_fs = 0;
			ve += len, qe += (op != NS_CIGAR_N);
		}
	}
	f = &feat[j++];
	f->vs = vs, f->ve = ve, f->qs = qs, f->qe = qe, f->phase = phase, f->n_fs = n_fs, f->type = MP_FEAT_EXON;
	if (has_stop) {
		f = &feat[j++];
		f->vs = ve_mRNA - 3, f->ve = ve_mRNA, f->qs = qe, f->qe = qe, f->phase = 0, f->n_fs = 0, f->type = MP_FEAT_STOP;
	}

	ctg = &mi->nt->ctg[r->vid>>1];
	snprintf(buf, 16, "%.6ld", (long)id);
	vs = r->vid&1? ctg->len - ve_mRNA : r->vs;
	ve = r->vid&1? ctg->len - r->vs   : ve_mRNA;
	mp_sprintf_lite(s, "%s\tminiprot\tmRNA\t%d\t%d\t%d\t%c\t.\tID=%s%s", ctg->name, (int)vs + 1, (int)ve, r->p->dp_max, "+-"[r->vid&1], gff_prefix, buf);
	if (tot_fs > 0) mp_sprintf_lite(s, ";Frameshift=true");
	snprintf(dec, 16, "%.4f", (double)r->p->n_iden * 3 / r->p->clen);
	mp_sprintf_lite(s, ";Identity=%s", dec);
	snprintf(dec, 16, "%.4f", (double)r->p->n_plus * 3 / r->p->clen);
	mp_sprintf_lite(s, ";Positive=%s", dec);
	mp_sprintf_lite(s, ";Target=%s %d %d\n", seq->name, r->qs + 1, r->qe);

	for (j = 0; j < n_feat; ++j) {
		f = &feat[j];
		vs = r->vid&1? ctg->len - f->ve : f->vs;
		ve = r->vid&1? ctg->len - f->vs : f->ve;
		mp_sprintf_lite(s, "%s\tminiprot\t%s\t%d\t%d\t0\t%c\t%d\tParent=%s%s", ctg->name, f->type == MP_FEAT_STOP? "stop_codon" : "CDS",
			(int)vs + 1, (int)ve, "+-"[r->vid&1], f->phase, gff_prefix, buf);
		if (f->n_fs > 0) mp_sprintf_lite(s, ";Frameshift=true");
		if (f->type == MP_FEAT_EXON)
			mp_sprintf_lite(s, ";Target=%s %d %d", seq->name, f->qs + 1, f->qe);
		mp_sprintf_lite(s, "\n");
	}
	kfree(km, feat);
}

void mp_write_output(kstring_t *s, void *km, const mp_idx_t *mi, const mp_bseq1_t *seq, const mp_reg1_t *r, const mp_mapopt_t *opt, int64_t id)
{
	s->l = 0;
	if (r == 0) {
		if (opt->flag&MP_F_SHOW_UNMAP)
			mp_write_paf(s, mi, seq, 0, opt->flag&MP_F_GFF);
	} else {
		mp_write_paf(s, mi, seq, r, opt->flag&MP_F_GFF);
		if (opt->flag&MP_F_GFF)
			mp_write_gff(s, km, mi, seq, r, opt->gff_prefix, id, seq->name);
	}
}
