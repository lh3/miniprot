#include <stdarg.h>
#include <string.h>
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

void mp_write_paf(kstring_t *s, const mp_idx_t *mi, const mp_bseq1_t *seq, const mp_reg1_t *r)
{
	const mp_ctg_t *ctg = &mi->nt->ctg[r->vid>>1];
	mp_sprintf_lite(s, "%s\t%d\t%d\t%d\t%c\t%s\t%d\t", seq->name, seq->l_seq, r->qs, r->qe, "+-"[r->vid&1], ctg->name, (int)ctg->len);
	if (r->vid&1) mp_sprintf_lite(s, "%d\t%d\t", (int)(ctg->len - r->ve), (int)(ctg->len - r->vs)); // FIXME: make it work for 64-bit integers
	else mp_sprintf_lite(s, "%d\t%d\t", (int)r->vs, (int)r->ve);
	if (r->p) {
		int32_t k;
		mp_sprintf_lite(s, "%d\t%d\tAS:i:%d\tnp:i:%d\tgl:i:%d\t", r->p->n_iden * 3, r->p->clen, r->p->dp_max, r->p->n_plus, r->p->glen);
		mp_sprintf_lite(s, "cg:Z:");
		for (k = 0; k < r->p->n_cigar; ++k)
			mp_sprintf_lite(s, "%d%c", r->p->cigar[k]>>4, NS_CIGAR_STR[r->p->cigar[k]&0xf]);
	} else mp_sprintf_lite(s, "%d\t%d", r->chn_sc, r->cnt);
	mp_sprintf_lite(s, "\n");
}
