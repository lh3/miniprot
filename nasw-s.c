#include <string.h>
#include "nasw.h"
#include "kalloc.h"
#include "mppriv.h"

#define NS_NEG_INF (-0x40000000)
#define ns_max(x, y) ((x) >= (y)? (x) : (y))

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

static void ns_dps_backtrack(void *km, const uint8_t *bk, int32_t nal, int32_t aal, uint32_t **cigar_, int32_t *n_cigar, int32_t *m_cigar)
{
	int32_t i = nal - 1, j = aal - 1, last = 0;
	uint32_t *cigar = *cigar_;
	while (i >= 0 && j >= 0) {
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
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_N, 1), --i;
			if (!ext) --j;
		} else if (state == 5) {
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_N, 1), --i;
			if (!ext) --j;
		} else if (state == 6) {
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_F, 1), --i;
		} else if (state == 7) {
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_G, 1), i -= 2;
		}
		last = state >= 1 && state <= 5 && ext? state : 0;
	}
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
void ns_dps_align_splice(void *km, const char *ns, int32_t nl, const char *as, int32_t al, int32_t asize, const int8_t *mat, int32_t q, int32_t e, int32_t r, int32_t f, int32_t cp)
{
	int32_t nal, aal, i, j, *G, *H, *I, *D;
	uint8_t *nas, *aas, *bk = 0;
	int8_t *nap, *acceptor, *donor;

	{ // generate nas[] and aas[]
		int32_t l;
		uint8_t codon;
		nal = nl + 1, aal = al;
		nas = Kmalloc(km, uint8_t, nal + aal);
		aas = nas + nal;
		for (j = 0; j < aal; ++j)
			aas[j] = mp_tab_aa20[(uint8_t)as[j]];
		memset(nas, 21, nal);
		for (i = l = 0, codon = 0; i < nl; ++i) {
			uint8_t c = mp_tab_nt4[(uint8_t)ns[i]];
			if (c < 4) {
				codon = (codon << 2 | c) & 0x3f;
				if (++l >= 3)
					nas[i+1] = mp_tab_codon[codon];
			} else codon = 0, l = 0;
		}
	}

	{ // generate nap[], donor[] and acceptor[]
		int32_t k, c;
		nap = Kmalloc(km, int8_t, nal * (asize + 2));
		donor = nap + nal * asize, acceptor = donor + nal;
		for (c = 0, k = 0; c < asize; ++c) {
			const int8_t *p = &mat[asize * c];
			for (i = 0; i < nal; ++i)
				nap[k++] = p[nas[i]];
		}
		for (i = 0; i < nal; ++i)
			donor[i] = acceptor[i] = -cp;
		for (i = 0; i < nl - 3; ++i) {
			int32_t t = 0;
			if (ns[i+1] == 2 && ns[i+2] == 3) t = 1;
			if (t && i + 3 < nl && (ns[i+3] == 0 || ns[i+3] == 2)) t = 2;
			donor[i+1] = t == 2? 0 : t == 1? -cp/2 : -cp;
		}
		for (i = 1; i < nl; ++i) {
			int32_t t = 0;
			if (ns[i-1] == 0 && ns[i] == 2) t = 1;
			if (t && i > 0 && (ns[i-2] == 1 || ns[i-2] == 3)) t = 2;
			acceptor[i+1] = t == 2? 0 : t == 1? -cp/2 : -cp;
		}
	}

	/*
	 * I(i,j) = max{ H(i,j-1) - q, I(i,j-1) } - e
	 * D(i,j) = max{ H(i-3,j) - q, D(i-3,j) } - e
	 * A(i,j) = max{ H(i-1,j)   - r - d(i-1), A(i-1,j) }
	 * B(i,j) = max{ H(i-1,j-1) - r - d(i),   B(i-1,j) }
	 * C(i,j) = max{ H(i-1,j-1) - r - d(i+1), C(i-1,j) }
	 * H(i,j) = max{ H(i-3,j-1) + s(i,j), I(i,j), D(i,j), H(i-2,j)-f, H(i-1,j)-f, A(i,j)-a(i), B(i,j)-a(i-2), C(i,j)-a(i-1) }
	 */
	{ // initialization
		int32_t A;
		H = Kmalloc(km, int32_t, nal * 4);
		G = H + nal, I = G + nal, D = I + nal;
		for (i = 0; i < nal * 4; ++i) H[i] = NS_NEG_INF;
		G[0] = G[1] = G[2] = 0;
		D[0] = D[1] = D[2] = -q;
		for (i = 3, A = -r; i < nal; ++i) {
			D[i] = D[i-3] - e;
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

				z |= G[i] - q >= I[i]? 0 : 1<<3;
				I[i] = ns_max(G[i] - q, I[i]) - e;
				y = h >= I[i]? y : 1;
				h = ns_max(h, I[i]);

				z |= H[i-3] - q >= D[i-3]? 0 : 1<<4;
				D[i] = ns_max(H[i-3] - q, D[i-3]) - e;
				y = h >= D[i]? y : 2;
				h = ns_max(h, D[i]);

				tmp = H[i-1] - r - donor[i-1];
				z |= tmp >= A? 0 : 1<<5;
				A = ns_max(tmp, A);
				y = h >= A - acceptor[i]? y : 3;
				h = ns_max(h, A - acceptor[i]);

				tmp = G[i-2] - r - donor[i];
				z |= tmp >= B? 0 : 1<<6;
				B = ns_max(tmp, B);
				y = h >= B - acceptor[i-2]? y : 4;
				h = ns_max(h, B - acceptor[i-2]);

				tmp = G[i-1] - r - donor[i+1];
				z |= tmp >= C? 0 : 1<<7;
				C = ns_max(tmp, C);
				y = h >= C - acceptor[i-1]? y : 5;
				h = ns_max(h, C - acceptor[i-1]);

				y = h >= H[i-2] - f? y : 6;
				h = ns_max(h, H[i-2] - f);

				y = h >= H[i-1] - f? y : 7;
				h = ns_max(h, H[i-1] - f);

				H[i] = h, bkj[i] = z | y;
			}
			swap = G, G = H, H = swap;
		}
	}

	// free
	kfree(km, bk);
	kfree(km, H);   // along with G[], D[] and I[]
	kfree(km, nap); // along with donor[] and acceptor[]
	kfree(km, nas); // along with aas[]
}
