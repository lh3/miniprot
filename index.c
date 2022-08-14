#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "mppriv.h"
#include "kalloc.h"
#include "kvec-km.h"
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
			d->ctg = Krealloc(0, mp_ctg_t, d->ctg, d->m_ctg);
		}
		c = &d->ctg[d->n_ctg++];
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

uint32_t *mp_idx_boff(const mp_ntdb_t *db, int32_t bbit)
{
	int32_t i;
	int64_t boff = 0;
	uint32_t *bo;
	bo = Kmalloc(0, uint32_t, db->n_ctg * 2);
	for (i = 0; i < db->n_ctg; ++i) {
		bo[i<<1|0] = boff;
		boff += (db->ctg[i].len + (1<<bbit) - 1) >> bbit;
		bo[i<<1|1] = boff;
		boff += (db->ctg[i].len + (1<<bbit) - 1) >> bbit;
	}
	assert(boff < UINT32_MAX);
	return bo;
}

int64_t mp_ntseq_get(const mp_ntdb_t *db, int32_t cid, int64_t st, int64_t en, int32_t rev, uint8_t *seq)
{
	int64_t i, s, e, k;
	if (cid >= db->n_ctg || cid < 0) return -1;
	if (en < 0 || en > db->ctg[cid].len) en = db->ctg[cid].len;
	s = db->ctg[cid].off + st;
	e = db->ctg[cid].off + en;
	if (!rev) {
		for (i = s, k = 0; i < e; ++i)
			seq[k++] = db->seq[i>>1] >> ((i&1) * 4) & 0xf;
	} else {
		for (i = e - 1, k = 0; i >= s; --i) {
			uint8_t c = db->seq[i>>1] >> ((i&1) * 4) & 0xf;
			seq[k++] = c >= 4? c : 3 - c;
		}
	}
	return k;
}

void mp_idx_proc_orf(void *km, const uint8_t *seq, int64_t st, int64_t en, int32_t kmer, int32_t smer, int32_t bbit, int64_t boff, mp64_v *a)
{
	int64_t i;
	int32_t l;
	uint32_t x, mask_k = (1U<<kmer*4) - 1, mask_s = (1U<<smer*4) - 1;
	for (i = st, l = 0, x = 0; i < en; i += 3) {
		uint8_t codon = seq[i]<<4 | seq[i+1]<<2 | seq[i+2];
		x = (x<<4 | mp_tab_codon13[codon]) & mask_k;
		if (++l >= kmer) {
			int32_t sel = 0;
			uint32_t y = mp_hash32_mask(x, mask_k);
			if (kmer > smer) {
				int32_t j;
				uint32_t m = UINT32_MAX;
				for (j = 0; j < (kmer - smer) << 2; j += 2) {
					uint32_t z = y>>j & mask_s;
					m = m < z? m : z;
				}
				sel = (m == (y>>(kmer-smer)*2 & mask_s));
			} else sel = 1;
			if (sel)
				kv_push(uint64_t, km, *a, (uint64_t)y<<32 | ((st>>bbit) + boff));
		}
	}
}

void mp_idx_proc_seq(void *km, int64_t len, const uint8_t *seq, int32_t min_aa_len, int32_t kmer, int32_t smer, int32_t bbit, int64_t boff, mp64_v *a)
{
	uint8_t codon[3], p;
	int64_t i, j, e[3], k[3], l[3];
	a->n = 0;
	for (i = 0; i < 3; ++i)
		e[i] = -1, k[i] = l[i] = 0, codon[i] = 0;
	for (i = 0, p = 0; i < len; ++i, ++p) {
		if (p == 3) p = 0;
		if (seq[i] < 4) {
			codon[p] = (codon[p] << 2 | seq[i]) & 0x3f;
			if (++l[p] >= 3) {
				uint8_t aa = mp_tab_codon[(uint8_t)codon[p]];
				if (aa >= 20) {
					if (k[p] >= min_aa_len)
						mp_idx_proc_orf(km, seq, e[p] + 1 - k[p] * 3, e[p] + 1, kmer, smer, bbit, boff, a);
					k[p] = l[p] = 0, e[p] = -1;
				} else e[p] = i, ++k[p];
			}
		} else {
			if (k[p] >= min_aa_len)
				mp_idx_proc_orf(km, seq, e[p] + 1 - k[p] * 3, e[p] + 1, kmer, smer, bbit, boff, a);
			k[p] = l[p] = 0, e[p] = -1;
		}
	}
	for (i = 0; i < 3; ++i)
		if (k[p] >= min_aa_len)
			mp_idx_proc_orf(km, seq, e[p] + 1 - k[p] * 3, e[p] + 1, kmer, smer, bbit, boff, a);
	radix_sort_mp64(a->a, a->a + a->n);
	for (i = 1, j = 0; i < a->n; ++i)
		if (a->a[j] != a->a[i])
			a->a[++j] = a->a[i];
	a->n = ++j;
}

mp_idx_t *mp_index(const mp_idxopt_t *io, const char *fn)
{
	mp_ntdb_t *nt;
	mp_idx_t *mi;
	nt = mp_ntseq_read(fn);
	if (nt == 0) return 0;
	if (mp_verbose >= 3)
		fprintf(stderr, "[M::%s@%.3f] read %ld bases in %d contigs\n", __func__, mp_realtime(), (long)nt->l_seq, nt->n_ctg);
	mi = Kcalloc(0, mp_idx_t, 1);
	mi->nt = nt;
	mi->boff = mp_idx_boff(mi->nt, io->bbit);
	return mi;
}
