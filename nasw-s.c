#include <string.h>
#include <stdio.h>
#include "nasw.h"
#include "kalloc.h"

#define NS_NEG_INF (-0x40000000)
#define ns_max(x, y) ((x) >= (y)? (x) : (y))

void ns_rst_init(ns_rst_t *r)
{
	memset(r, 0, sizeof(*r));
}

static inline uint32_t *ns_push_cigar(void *km, int32_t *n_cigar, int32_t *m_cigar, uint32_t *cigar, uint32_t op, int32_t len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			(*m_cigar) += ((*m_cigar)>>1) + 8;
			cigar = Krealloc(km, uint32_t, cigar, *m_cigar);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

static void ns_s1_backtrack(void *km, const uint8_t *bk, int32_t nal, int32_t aal, uint32_t **cigar_, int32_t *n_cigar, int32_t *m_cigar)
{
	int32_t i = nal - 1, j = aal - 1, last = 0;
	uint32_t *cigar = *cigar_, tmp;
	while (i >= 1 && j >= 0) {
		int32_t x = bk[j * nal + i], state, ext;
		state = last == 0? x&7 : last;
		ext = state >= 1 && state <= 5? x>>(state+2)&1 : 0;
		if (state == 0) {
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_M, 1), i -= 3, --j;
		} else if (state == 1) {
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_I, 1), --j;
		} else if (state == 2) {
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_D, 1), i -= 3;
		} else if (state == 3) {
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_N, 1), --i;
		} else if (state == 4) {
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_U, 1), --i;
			if (!ext) --j;
		} else if (state == 5) {
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_V, 1), --i;
			if (!ext) --j;
		} else if (state == 6) {
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_F, 1), --i;
		} else if (state == 7) {
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_G, 2), i -= 2;
		}
		last = state >= 1 && state <= 5 && ext? state : 0;
	}
	if (j > 0) ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_I, j); // TODO: is this correct?
	if (i >= 1) {
		int32_t l = (i-1)/3*3, t = i - l;
		if (l > 0) ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_D, l); // TODO: is this correct?
		if (t == 2) ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_G, 2);
		else if (t == 1) ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_F, 1);
	}
	for (i = 0; i < (*n_cigar)>>1; ++i) // reverse CIGAR
		tmp = cigar[i], cigar[i] = cigar[(*n_cigar) - 1 - i], cigar[(*n_cigar) - 1 - i] = tmp;
	*cigar_ = cigar;
}

/*
 * M(i,j) = max{ M(i-3,j-1), I(i-3,j-1), D(i-3,j-1), A(i-1,j-1)-a(i-1), B(i-3,j-1)-a(i-3), C(i-2,j-1)-a(i-2) } + s(i,j)
 *
 * H(-1,-1) = H(0,-1) = H(1,-1) = 0; for i<-1, H(i,-1) = -inf
 * I(-1,0) = I(0,0) = I(1,0) = -q - e
 * D(-1,-1) = D(0,-1) = D(1,-1) = -q
 * A(0,-1) = -r - d
 * B(1,-1) = C(2,-1) = -inf
 *
 * H(i,j)   = max{ H(i-3,j-1) + s(i,j), I(i,j), D(i,j), H(i-2,j)-f, H(i-1,j)-f, A(i,j)-a(i), B(i,j)-a(i-2), C(i,j)-a(i-1) }
 * I(i,j+1) = max{ H(i,j) - q, I(i,j) } - e
 * D(i+3,j) = max{ H(i,j) - q, D(i,j) } - e
 * A(i+1,j) = max{ H(i,j)   - r - d(i),   A(i,j) }
 * B(i+1,j) = max{ H(i,j-1) - r - d(i+1), B(i,j) }
 * C(i+1,j) = max{ H(i,j-1) - r - d(i+2), C(i,j) }
 */
void ns_splice_s1(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, ns_rst_t *r)
{
	int32_t nal, aal, i, j, *mem_H, *G, *H, *I, *D;
	uint8_t *nas, *aas, *bk = 0;
	int8_t *nap, *acceptor, *donor;

	r->n_cigar = 0;

	{ // generate nas[], aas[], donor[], acceptor[] and nap[]
		int32_t l, k, c;
		uint8_t codon;
		nal = nl + 1, aal = al;
		nas = Kmalloc(km, uint8_t, nal + aal + (nal + 1) * 2 + nal * opt->asize);
		aas = nas + nal;
		donor = (int8_t*)aas + aal, acceptor = donor + nal + 1;
		nap = acceptor + nal + 1;
		for (j = 0; j < aal; ++j) // generate aas[]
			aas[j] = opt->aa20[(uint8_t)as[j]];
		for (i = 0; i < nl; ++i) // nt4 encoding of ns[] for computing donor[] and acceptor[]
			nas[i] = opt->nt4[(uint8_t)ns[i]];
		for (i = 0; i < nal + 1; ++i)
			donor[i] = acceptor[i] = opt->nc;
		for (i = 0; i < nl - 3; ++i) { // generate donor[]
			int32_t t = 0;
			if (nas[i+1] == 2 && nas[i+2] == 3) t = 1;
			if (t && i + 3 < nl && (nas[i+3] == 0 || nas[i+3] == 2)) t = 2;
			donor[i+1] = t == 2? 0 : t == 1? opt->nc/2 : opt->nc;
		}
		for (i = 1; i < nl; ++i) { // generate acceptor[]
			int32_t t = 0;
			if (nas[i-1] == 0 && nas[i] == 2) t = 1;
			if (t && i > 0 && (nas[i-2] == 1 || nas[i-2] == 3)) t = 2;
			acceptor[i+1] = t == 2? 0 : t == 1? opt->nc/2 : opt->nc;
		}
		memset(nas, opt->aa20['X'], nal);
		for (i = l = 0, codon = 0; i < nl; ++i) { // generate the real nas[]
			uint8_t c = opt->nt4[(uint8_t)ns[i]];
			if (c < 4) {
				codon = (codon << 2 | c) & 0x3f;
				if (++l >= 3)
					nas[i+1] = opt->codon[codon];
			} else codon = 0, l = 0;
		}
		nap = Kmalloc(km, int8_t, nal * opt->asize);
		for (c = 0, k = 0; c < opt->asize; ++c) { // generate nap[]
			const int8_t *p = &opt->sc[opt->asize * c];
			for (i = 0; i < nal; ++i)
				nap[k++] = p[nas[i]];
		}
		//for (i = 0; i < nal; ++i) putchar(ns_tab_aa_i2c[nas[i]]); putchar('\n');
	}

	/*
	 * I(i,j) = max{ H(i,j-1) - q, I(i,j-1) } - e
	 * D(i,j) = max{ H(i-3,j) - q, D(i-3,j) } - e
	 * A(i,j) = max{ H(i-1,j)   - r - d(i-1), A(i-1,j) }
	 * B(i,j) = max{ H(i-1,j-1) - r - d(i),   B(i-1,j) }
	 * C(i,j) = max{ H(i-1,j-1) - r - d(i+1), C(i-1,j) }
	 * H(i,j) = max{ H(i-3,j-1) + s(i,j), I(i,j), D(i,j), H(i-1,j)-f, H(i-2,j)-f, A(i,j)-a(i), B(i,j)-a(i-2), C(i,j)-a(i-1) }
	 */
	{ // initialization. FIXME: I[] is not initialized correctly
		int32_t A;
		mem_H = H = Kmalloc(km, int32_t, nal * 4);
		G = H + nal, I = G + nal, D = I + nal;
		for (i = 0; i < nal * 4; ++i) H[i] = NS_NEG_INF;
		G[0] = 0, G[1] = G[2] = -opt->fs;
		D[0] = -opt->go, D[1] = D[2] = -opt->fs - opt->go;
		for (i = 3, A = -opt->io; i < nal; ++i) {
			D[i] = D[i-3] - opt->ge;
			G[i] = ns_max(D[i], A - acceptor[i]);
		}
	}
	{ // core loop
		bk = Kmalloc(km, uint8_t, nal * aal);
		for (j = 0; j < aal; ++j) {
			uint8_t *bkj = &bk[j * nal];
			int8_t *ms = &nap[aas[j] * nal];
			int32_t A, B, C, *swap;
			A = B = C = D[0] = D[1] = D[2] = NS_NEG_INF; // FIXME: this is not correct
			for (i = 3; i < nal; ++i) {
				uint8_t z = 0, y = 0;
				int32_t h = G[i-3] + ms[i], tmp;
				//printf("%d\t%d\n", i, donor[i]);

				z |= G[i] - opt->go >= I[i]? 0 : 1<<3;
				I[i] = ns_max(G[i] - opt->go, I[i]) - opt->ge;
				y = h >= I[i]? y : 1;
				h = ns_max(h, I[i]);

				z |= H[i-3] - opt->go >= D[i-3]? 0 : 1<<4;
				D[i] = ns_max(H[i-3] - opt->go, D[i-3]) - opt->ge;
				y = h >= D[i]? y : 2;
				h = ns_max(h, D[i]);

				tmp = H[i-1] - opt->io - donor[i-1];
				z |= tmp >= A? 0 : 1<<5;
				A = ns_max(tmp, A);
				y = h >= A - acceptor[i]? y : 3;
				h = ns_max(h, A - acceptor[i]);

				tmp = G[i-1] - opt->io - donor[i];
				z |= tmp >= B? 0 : 1<<6;
				B = ns_max(tmp, B);
				y = h >= B - acceptor[i-2]? y : 4;
				h = ns_max(h, B - acceptor[i-2]);

				tmp = G[i-1] - opt->io - donor[i+1];
				z |= tmp >= C? 0 : 1<<7;
				C = ns_max(tmp, C);
				y = h >= C - acceptor[i-1]? y : 5;
				h = ns_max(h, C - acceptor[i-1]);

				y = h >= H[i-1] - opt->fs? y : 6;
				h = ns_max(h, H[i-1] - opt->fs);

				y = h >= H[i-2] - opt->fs? y : 7;
				h = ns_max(h, H[i-2] - opt->fs);

				H[i] = h, bkj[i] = z | y;
			}
			swap = G, G = H, H = swap;
		}
		r->score = G[nal - 1];
	}

	// free
	kfree(km, mem_H); // H[], G[], D[] and I[]
	kfree(km, nas);   // along with aas[]

	ns_s1_backtrack(km, bk, nal, aal, &r->cigar, &r->n_cigar, &r->m_cigar);
	kfree(km, bk);
}
