#include "mppriv.h"
#include "kalloc.h"

#define MP_NEG_INF (-0x40000000)
#define mp_max(x, y) ((x) >= (y)? (x) : (y))

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
 * A(i+1,j) = max{ H(i,j)     - r - d(i), A(i,j) }
 * B(i+1,j) = max{ H(i-1,j-1) - r - d(i), B(i,j) }
 * C(i+1,j) = max{ H(i-2,j-1) - r - d(i), C(i,j) }
 */
void mp_dps_align_splice(void *km, const char *ns, int32_t nl, const char *as, int32_t al, int32_t asize, const int8_t *mat, int32_t q, int32_t e, int32_t r, int32_t f, int32_t cp)
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
 * B(i,j) = max{ H(i-2,j-1) - r - d(i-1), B(i-1,j) }
 * C(i,j) = max{ H(i-3,j-1) - r - d(i-1), C(i-1,j) }
 * H(i,j) = max{ H(i-3,j-1) + s(i,j), I(i,j), D(i,j), H(i-2,j)-f, H(i-1,j)-f, A(i,j)-a(i), B(i,j)-a(i-2), C(i,j)-a(i-1) }
 */
	{ // initialization
		int32_t A;
		H = Kmalloc(km, int32_t, nal * 4);
		G = H + nal, I = G + nal, D = I + nal;
		for (i = 0; i < nal * 4; ++i) H[i] = MP_NEG_INF;
		G[0] = G[1] = G[2] = 0;
		D[0] = D[1] = D[2] = -q;
		for (i = 3, A = -r; i < nal; ++i) {
			D[i] = D[i-3] - e;
			G[i] = mp_max(D[i], A - acceptor[i]);
		}
	}
	{ // core loop
		for (j = 0; j < aal; ++j) {
			int8_t *ms = &nap[aas[j] * nal];
			int32_t A, B, C, *swap;
			A = B = C = D[0] = D[1] = D[2] = MP_NEG_INF; // FIXME: this is not correct
			for (i = 3; i < nal; ++i) {
				int32_t M, rd = r - donor[i], gp, sl, fs;
				I[i] = mp_max(G[i] - q, I[i]) - e;
				D[i] = mp_max(H[i-3] - q, D[i-3]) - e;
				A = mp_max(H[i-1] - rd, A);
				B = mp_max(G[i-2] - rd, B);
				C = mp_max(G[i-3] - rd, C);
				M = G[i-3] + ms[i];
				gp = mp_max(I[i], D[i]);
				sl = mp_max(A - acceptor[i], B - acceptor[i-2]);
				sl = mp_max(sl, C - acceptor[i-1]);
				fs = mp_max(H[i-2], H[i-1]) - f;
				M = mp_max(M, gp);
				M = mp_max(M, sl);
				M = mp_max(M, fs);
				H[i] = M;
			}
			swap = G, G = H, H = swap;
		}
	}

	// free
	kfree(km, H);   // along with G[], D[] and I[]
	kfree(km, nap); // along with donor[] and acceptor[]
	kfree(km, nas); // along with aas[]
}
