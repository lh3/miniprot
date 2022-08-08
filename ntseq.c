#include <stdio.h>
#include <zlib.h>
#include "mppriv.h"
#include "kalloc.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

mp_ntdb_t *mp_ntseq_read(const char *fn)
{
	gzFile fp;
	kseq_t *ks;
	mp_ntdb_t *d = 0;
	int64_t off = 0;

	fp = gzopen(fn, "r");
	if (fp == 0) return 0;
	ks = kseq_init(fp);

	d = Kcalloc(0, mp_ntdb_t, 1);
	while (kseq_read(ks) >= 0) {
		int64_t i, ltmp;
		mp_ctg_t *c;

		// update mp_ntdb_t::ctg
		if (d->n_ctg == d->m_ctg) {
			d->m_ctg += (d->m_ctg>>1) + 16;
			d->ctg = Krealloc(0, mp_ctg_t, d->m_ctg);
		}
		c = d->ctg[d->n_ctg++];
		c->name = mp_strdup(ks->name.s);
		c->off = off;
		c->len = ks->seq.l;

		// update mp_ntdb_t::seq
		ltmp = (d->l_seq + ks->seq.l + 1) >> 1 << 1;
		if (ltmp > d->m_seq) {
			int64_t oldm = d->m_seq;
			d->m_seq = ltmp;
			kroundup64(d->m_seq);
			d->seq = Krealloc(0, uint8_t, d->seq, d->m_seq);
			memset(&d->seq[oldm>>1], 0, (d->m_seq - oldm) >> 1);
		}
		for (i = 0; i < ks->seq.l; ++i, ++off) {
			uint8_t b = mp_tab_nt4[(uint8_t)ks->seq.s[i]];
			d->seq[off >> 1] |= b << (off&1) * 4;
		}
	}

	kseq_destroy(ks);
	gzclose(fp);
	return d;
}

int64_t mp_ntseq_get(const mp_ntdb_t *db, int32_t cid, int64_t st, int64_t en, int32_t rev, uint8_t *seq)
{
}
