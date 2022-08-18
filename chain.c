#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "mppriv.h"
#include "kalloc.h"

static int64_t mp_chain_bk_end(int32_t max_drop, const mp128_t *z, const int32_t *f, const int64_t *p, int32_t *t, int64_t k)
{
	int64_t i = z[k].y, end_i = -1, max_i = i;
	int32_t max_s = 0;
	if (i < 0 || t[i] != 0) return i;
	do {
		int32_t s;
		t[i] = 2;
		end_i = i = p[i];
		s = i < 0? z[k].x : (int32_t)z[k].x - f[i];
		if (s > max_s) max_s = s, max_i = i;
		else if (max_s - s > max_drop) break;
	} while (i >= 0 && t[i] == 0);
	for (i = z[k].y; i >= 0 && i != end_i; i = p[i]) // reset modified t[]
		t[i] = 0;
	return max_i;
}

static uint64_t *mp_chain_backtrack(void *km, int64_t n, const int32_t *f, const int64_t *p, int32_t *v, int32_t *t, int32_t min_cnt, int32_t min_sc, int32_t max_drop, int32_t *n_u_, int32_t *n_v_)
{
	mp128_t *z;
	uint64_t *u;
	int64_t i, k, n_z, n_v;
	int32_t n_u;

	*n_u_ = *n_v_ = 0;
	for (i = 0, n_z = 0; i < n; ++i) // precompute n_z
		if (f[i] >= min_sc) ++n_z;
	if (n_z == 0) return 0;
	KMALLOC(km, z, n_z);
	for (i = 0, k = 0; i < n; ++i) // populate z[]
		if (f[i] >= min_sc) z[k].x = f[i], z[k++].y = i;
	radix_sort_mp128x(z, z + n_z);

	memset(t, 0, n * 4);
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // precompute n_u
		if (t[z[k].y] == 0) {
			int64_t n_v0 = n_v, end_i;
			int32_t sc;
			end_i = mp_chain_bk_end(max_drop, z, f, p, t, k);
			for (i = z[k].y; i != end_i; i = p[i])
				++n_v, t[i] = 1;
			sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
			if (sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt)
				++n_u;
			else n_v = n_v0;
		}
	}
	KMALLOC(km, u, n_u);
	memset(t, 0, n * 4);
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // populate u[]
		if (t[z[k].y] == 0) {
			int64_t n_v0 = n_v, end_i;
			int32_t sc;
			end_i = mp_chain_bk_end(max_drop, z, f, p, t, k);
			for (i = z[k].y; i != end_i; i = p[i])
				v[n_v++] = i, t[i] = 1;
			sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
			if (sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt)
				u[n_u++] = (uint64_t)sc << 32 | (n_v - n_v0);
			else n_v = n_v0;
		}
	}
	kfree(km, z);
	assert(n_v < INT32_MAX);
	*n_u_ = n_u, *n_v_ = n_v;
	return u;
}

static uint64_t *compact_a(void *km, int32_t n_u, uint64_t *u, int32_t n_v, int32_t *v, uint64_t *a)
{
	mp128_t *w;
	uint64_t *b, *u2;
	int64_t i, j, k;

	// write the result to b[]
	KMALLOC(km, b, n_v);
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			b[k++] = a[v[k0 + (ni - j - 1)]];
	}
	kfree(km, v);

	// sort u[] and a[] by the target position, such that adjacent chains may be joined
	KMALLOC(km, w, n_u);
	for (i = k = 0; i < n_u; ++i) {
		w[i].x = b[k]>>32, w[i].y = (uint64_t)k<<32|i;
		k += (int32_t)u[i];
	}
	radix_sort_mp128x(w, w + n_u);
	KMALLOC(km, u2, n_u);
	for (i = k = 0; i < n_u; ++i) {
		int32_t j = (int32_t)w[i].y, n = (int32_t)u[j];
		u2[i] = u[j];
		memcpy(&a[k], &b[w[i].y>>32], n * sizeof(*a));
		k += n;
	}
	memcpy(u, u2, n_u * 8);
	memcpy(b, a, k * sizeof(*a)); // write _a_ to _b_ and deallocate _a_ because _a_ is oversized, sometimes a lot
	kfree(km, a); kfree(km, w); kfree(km, u2);
	return b;
}

static inline int32_t comput_sc(uint64_t ai, uint64_t aj, int32_t max_dist_x, int32_t max_dist_y, int32_t bw, int32_t is_spliced, int32_t bbit, int32_t kmer)
{
	int32_t dq = (int32_t)ai - (int32_t)aj, dq3 = dq * 3, dr3, dd, sc;
	if (dq <= 0 || dq3 > max_dist_x) return INT32_MIN;
	if (dq > max_dist_y) return INT32_MIN;
	if (bbit > 0) {
		int32_t bs = 1<<bbit;
		dr3 = ((ai>>32) - (aj>>32)) << bbit;
		if (dq3 >= dr3 - bs && dq3 <= dr3 + bs) dd = 0;
		else if (dq3 < dr3 - bs) dd = dr3 - bs - dq3;
		else dd = dq3 - (dr3 + bs);
	} else {
		dr3 = (ai>>32) - (aj>>32);
		dd = dr3 > dq3? dr3 - dq3 : dq3 - dr3;
	}
//	if (ai>>32 == 25318 && (uint32_t)ai == 96 && aj>>32 == 25306 && (uint32_t)aj == 87) printf("here: %d,%d\n", dd, bw);
	if (dd > bw) return INT32_MIN; // dd is the min possible gap size
	if (bbit > 0) {
		sc = kmer < dq? kmer : dq;
	} else if (kmer <= dq && kmer * 3 <= dr3) {
		sc = kmer;
	} else {
		int32_t dr = dr3 / 3, q = dr3 - dr * 3;
		int32_t dg = dr < dq? dr : dq;
		sc = dg < kmer? dg : kmer;
		if (q != 0) --sc; // frameshift
	}
	if (dd > 0) {
		float lin_pen, log_pen;
		lin_pen = (float)dd;
		log_pen = dd >= 1? mp_log2(dd + 1) : 0.0f; // mp_log2() only works for dd>=2
		if (is_spliced) {
			if (dr3 > dq3) sc -= (int)(lin_pen < log_pen? lin_pen : log_pen);
			else sc -= (int)(lin_pen + .5f * log_pen);
		} else sc -= (int)(lin_pen + .5f * log_pen);
	}
	return sc;
}

/* Input:
 *   a: blockId << 32 | queryPos
 * Output:
 *   n_u: #chains
 *   u[]: score<<32 | #anchors (sum of lower 32 bits of u[] is the returned length of a[])
 * input a[] is deallocated on return
 */
uint64_t *mp_chain(int32_t max_dist_x, int32_t max_dist_y, int32_t bw, int32_t max_skip, int32_t max_iter, int32_t min_cnt, int32_t min_sc,
				   int32_t is_spliced, int32_t kmer, int32_t bbit, int64_t n, uint64_t *a, int32_t *n_u_, uint64_t **_u, void *km)
{ // TODO: make sure this works when n has more than 32 bits
	int32_t *f, *t, *v, n_u, n_v, mmax_f = 0, max_drop = bw;
	int64_t *p, i, j, st = 0, n_iter = 0;
	uint64_t *u;

	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) {
		kfree(km, a);
		return 0;
	}
	if (max_dist_x < bw) max_dist_x = bw;
	if (max_dist_y < bw && !is_spliced) max_dist_y = bw;
	if (is_spliced) max_drop = INT32_MAX;
	KMALLOC(km, p, n);
	KMALLOC(km, f, n);
	KMALLOC(km, v, n);
	KCALLOC(km, t, n);

	// fill the score and backtrack arrays
	for (i = 0; i < n; ++i) {
		int64_t max_j = -1;
		int32_t max_f = kmer, n_skip = 0;
		while (st < i && ((a[i]>>32) - (a[st]>>32)) << bbit > max_dist_x) ++st;
		if (i - st > max_iter) st = i - max_iter;
		for (j = i - 1; j >= st; --j) {
			int32_t sc;
			sc = comput_sc(a[i], a[j], max_dist_x, max_dist_y, bw, is_spliced, bbit, kmer);
			++n_iter;
			if (sc == INT32_MIN) continue;
			sc += f[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
				if (n_skip > 0) --n_skip;
			} else if (t[j] == (int32_t)i) {
				if (++n_skip > max_skip)
					break;
			}
			if (p[j] >= 0) t[p[j]] = i;
		}
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		if (mmax_f < max_f) mmax_f = max_f;
	}

	u = mp_chain_backtrack(km, n, f, p, v, t, min_cnt, min_sc, max_drop, &n_u, &n_v);
	*n_u_ = n_u, *_u = u; // NB: note that u[] may not be sorted by score here
	kfree(km, p); kfree(km, f); kfree(km, t);
	if (n_u == 0) {
		kfree(km, a); kfree(km, v);
		return 0;
	}
	return compact_a(km, n_u, u, n_v, v, a);
}
