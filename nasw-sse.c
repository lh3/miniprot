#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "nasw.h"
#include "kalloc.h"

#if defined(__SSE4_1__)
#include <smmintrin.h>
#elif defined(__SSE2__)
#include <xmmintrin.h>
#elif defined(__ARM_NEON)
#include "s2n-lite.h"
#endif

/*
 * I(i,j) = max{ H(i,j-1) - q, I(i,j-1) } - e
 * D(i,j) = max{ H(i-3,j) - q, D(i-3,j) } - e
 * A(i,j) = max{ H(i-1,j)   - r - d(i-1), A(i-1,j) }
 * B(i,j) = max{ H(i-1,j-1) - r - d(i),   B(i-1,j) }
 * C(i,j) = max{ H(i-1,j-1) - r - d(i+1), C(i-1,j) }
 * H(i,j) = max{ H(i-3,j-1) + s(i,j), I(i,j), D(i,j), H(i-1,j-1)-f, H(i-2,j-1)-f, H(i-1,j)-f, H(i-2,j)-f, A(i,j)-a(i), B(i,j)-a(i-2), C(i,j)-a(i-1) }
 */

static __m128i *ns_alloc16(void *km, size_t n, uint8_t **mem)
{
	*mem = Kmalloc(km, uint8_t, (16 * n + 31) / 16 * 16);
	return (__m128i*)(((size_t)(*mem) + 15) / 16 * 16);
}

static void ns_backtrack(void *km, int32_t vs, const __m128i *tb, int32_t nl, int32_t al, uint32_t **cigar_, int32_t *n_cigar, int32_t *m_cigar)
{
	int32_t i = nl - 1, j = al - 1, last = 0, slen = (al + vs - 1) / vs;
	uint32_t *cigar = *cigar_, tmp;
	assert(vs == 4 || vs == 8);
	while (i >= 2 && j >= 0) {
		const __m128i *tbi = &tb[i * slen];
		int32_t x = *((int16_t*)&tbi[j%slen] + j/slen);
		int32_t state, ext;
		if (vs == 4) x = *((int32_t*)&tbi[j%slen] + j/slen);
		if (x>>9&1) x = 1|x>>4<<4;
		state = last == 0? x&0xf : last;
		ext = state >= 1 && state <= 5? x>>(state+3)&1 : 0;
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
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_F, 2), i -= 2;
		} else if (state == 8) {
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_G, 1), --i, --j;
		} else if (state == 9) {
			cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_G, 2), i -= 2, --j;
		}
		last = state >= 1 && state <= 5 && ext? state : 0;
	}
	if (j > 0) cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_I, j);
	if (i >= 0) {
		int32_t l = (i+1) / 3 * 3, t = (i+1) % 3;
		if (l > 0) cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_D, l);
		if (t != 0) cigar = ns_push_cigar(km, n_cigar, m_cigar, cigar, NS_CIGAR_F, t);
	}
	for (i = 0; i < (*n_cigar)>>1; ++i) // reverse CIGAR
		tmp = cigar[i], cigar[i] = cigar[(*n_cigar) - 1 - i], cigar[(*n_cigar) - 1 - i] = tmp;
	//fprintf(stderr, "%d\t", nl); for (i = 0; i < (*n_cigar); ++i) fprintf(stderr, "%d%c", cigar[i]>>4, NS_CIGAR_STR[cigar[i]&0xf]); fputc('\n', stderr);
	*cigar_ = cigar;
}

static void ns_prep_nas(const char *ns, int32_t nl, const ns_opt_t *opt, uint8_t *nas)
{
	int32_t i, l;
	uint8_t codon;
	memset(nas, opt->aa20['X'], nl);
	for (i = l = 0, codon = 0; i < nl; ++i) { // generate the real nas[]
		uint8_t c = opt->nt4[(uint8_t)ns[i]];
		if (c < 4) {
			codon = (codon << 2 | c) & 0x3f;
			if (++l >= 3)
				nas[i] = opt->codon[codon];
		} else codon = 0, l = 0;
	}
}

static uint8_t *ns_prep_seq(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, const uint8_t *ss, uint8_t **aas_, int8_t **donor_, int8_t **acceptor_)
{
	int32_t i, j;
	uint8_t *nas, *aas;
	int8_t *donor, *acceptor;
	nas = Kmalloc(km, uint8_t, nl + al + (nl + 1) * 2); // nas and aas are allocated together
	*aas_ = aas = nas + nl;
	*donor_ = donor = (int8_t*)aas + al, *acceptor_ = acceptor = donor + nl + 1;
	for (j = 0; j < al; ++j) // generate aas[]
		aas[j] = opt->aa20[(uint8_t)as[j]];
	for (i = 0; i < nl; ++i) // nt4 encoding of ns[] for computing donor[] and acceptor[]
		nas[i] = opt->nt4[(uint8_t)ns[i]];
	for (i = 0; i < nl + 1; ++i)
		donor[i] = acceptor[i] = opt->sp[3];
	for (i = 0; i < nl - 3; ++i) { // generate donor[]
		int32_t t = 3;
		if (nas[i+1] == 2 && nas[i+2] == 3) // GT.
			t = i + 3 < nl && (nas[i+3] == 0 || nas[i+3] == 2)? (nas[i] == 2? -1 : 4) : 0;
		else if (nas[i+1] == 2 && nas[i+2] == 1 && nas[i] == 2) t = 1; // GC.
		else if (nas[i+1] == 0 && nas[i+2] == 3) t = 2; // AT.
		donor[i] = t < 0? 0 : opt->sp[t];
	}
	for (i = 1; i < nl; ++i) { // generate acceptor[]
		int32_t t = 3, penY = 0;
		if (nas[i-1] == 0 && nas[i] == 2) { // .AG
			t = i >= 2 && (nas[i-2] == 1 || nas[i-2] == 3)? -1 : 0;
			for (j = i - 4; j >= 0 && j > i - 7; --j)
				if (nas[j] != 1 && nas[j] != 3) penY += opt->sp[5];
		} else if (nas[i-1] == 0 && nas[i] == 1) t = 2; // .AC
		acceptor[i] = t < 0? 0 : opt->sp[t];
		if (t == -1 || t == 0) acceptor[i] += penY;
	}
	if (ss) { // NB: ss[] uses the offset rule to specify a splice site; donor/acceptor[] uses a different rule. The ss[] way is better but too tricky to change now
		for (i = 1; i < nl; ++i) {
			if (ss[i] == 0xff) { // score not set
				donor[i-1] -= opt->sp_null_bonus;
				acceptor[i-1] -= opt->sp_null_bonus;
			} else if (ss[i]&1) { // acceptor
				acceptor[i-1] -= (int8_t)(ss[i]>>1) - (int8_t)NS_SPSC_OFFSET;
			} else { // donor
				donor[i-1] -= (int8_t)(ss[i]>>1) - (int8_t)NS_SPSC_OFFSET;
			}
		}
	}
	ns_prep_nas(ns, nl, opt, nas);
	return nas;
}

static uint8_t *ns_prep_seq_left(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, const uint8_t *ss, uint8_t **aas_, int8_t **donor_, int8_t **acceptor_)
{
	int32_t i, j;
	uint8_t *nas, *aas, tmp;
	int8_t *donor, *acceptor;
	nas = Kmalloc(km, uint8_t, nl + al + (nl + 1) * 2); // nas and aas are allocated together
	*aas_ = aas = nas + nl;
	*donor_ = donor = (int8_t*)aas + al, *acceptor_ = acceptor = donor + nl + 1;
	for (j = 0; j < al; ++j) // generate aas[]. NB: this is reversed
		aas[al - 1 - j] = opt->aa20[(uint8_t)as[j]];
	for (i = 0; i < nl; ++i) // nt4 encoding of ns[] for computing donor[] and acceptor[]. NB: this is reversed
		nas[nl - 1 - i] = opt->nt4[(uint8_t)ns[i]];
	for (i = 0; i < nl + 1; ++i)
		donor[i] = acceptor[i] = opt->sp[3];
	for (i = 0; i < nl - 3; ++i) { // generate donor[] (actually acceptor as the string is reversed)
		int32_t t = 3, penY = 0;
		if (nas[i+1] == 2 && nas[i+2] == 0) { // GA. (the reverse of .AG)
			t = i + 3 < nl && (nas[i+3] == 1 || nas[i+3] == 3)? -1 : 0;
			for (j = i + 5; j < nl && j < i + 8; ++j)
				if (nas[j] != 1 && nas[j] != 3) penY += opt->sp[5];
		} else if (nas[i+1] == 1 && nas[i+2] == 0) t = 2; // CA.
		donor[i] = t < 0? 0 : opt->sp[t];
		if (t == -1 || t == 0) donor[i] += penY;
	}
	for (i = 1; i < nl; ++i) { // generate acceptor[] (actually donor)
		int32_t t = 3;
		if (nas[i-1] == 3 && nas[i] == 2) // .TG (the reverse of GT.)
			t = i >= 2 && (nas[i-2] == 0 || nas[i-2] == 2)? (i+1 < nl && nas[i+1] == 2? -1 : 4) : 0;
		else if (nas[i-1] == 1 && nas[i] == 2 && i+1 < nl && nas[i+1] == 1) t = 1; // .CG
		else if (nas[i-1] == 3 && nas[i] == 0) t = 2; // .TA
		acceptor[i] = t < 0? 0 : opt->sp[t];
	}
	if (ss) {
		for (i = 0; i < nl; ++i) {
			if (ss[i] == 0xff) {
				donor[nl - i - 1] -= opt->sp_null_bonus;
				acceptor[nl - i - 1] -= opt->sp_null_bonus;
			} else if (ss[i]&1) {
				donor[nl - i - 1] -= (int8_t)(ss[i]>>1) - (int8_t)NS_SPSC_OFFSET;
			} else {
				acceptor[nl - i - 1] -= (int8_t)(ss[i]>>1) - (int8_t)NS_SPSC_OFFSET;
			}
		}
	}
	ns_prep_nas(ns, nl, opt, nas);
	for (i = 0; i < nl>>1; ++i) // reverse
		tmp = nas[i], nas[i] = nas[nl - 1 - i], nas[nl - 1 - i] = tmp;
	memmove(nas + 2, nas, nl - 2);
	nas[0] = nas[1] = opt->aa20['X'];
	return nas;
}

#define ns_gen_prof(INT_TYPE, _km, aas, al, opt, neg_inf, _mem_ap, _ap) do { \
	INT_TYPE *t; \
	int32_t a, p = 16 / sizeof(INT_TYPE), slen = (al + p - 1) / p; \
	*(_ap) = ns_alloc16(_km, slen * opt->asize, _mem_ap); \
	t = (INT_TYPE*)*(_ap); \
	for (a = 0; a < opt->asize; ++a) { \
		int32_t i, k, nlen = slen * p; \
		const int8_t *ma = opt->sc + a * opt->asize;\
		for (i = 0; i < slen; ++i) \
			for (k = i; k < nlen; k += slen) \
				*t++ = (k >= al? neg_inf : ma[aas[k]]); \
	} \
} while (0)

#define sse_gen(func, suf) _mm_##func##_##suf

#define NS_GEN_VAR(_suf) \
	int32_t i, j, is_ext = !!(opt->flag&(NS_F_EXT_LEFT|NS_F_EXT_RIGHT)), slen = (al + vsize - 1) / vsize; /* segment length */ \
	uint8_t *nas, *aas, *mem_ap, *mem_H, *mem_tb = 0; \
	int8_t *donor, *acceptor; \
	__m128i *ap, *tb = 0, *H0, *H, *H1, *H2, *H3, *D, *D1, *D2, *D3, *A, *B, *C, *Hmax; \
	__m128i go, ge, goe, io, fs; \

#define NS_GEN_PREPARE(_suf) \
	r->n_cigar = 0, r->nt_len = nl, r->aa_len = al, r->score = INT32_MIN; \
	if (opt->flag & NS_F_EXT_LEFT) nas = ns_prep_seq_left(km, ns, nl, as, al, opt, ss, &aas, &donor, &acceptor); \
	else nas = ns_prep_seq(km, ns, nl, as, al, opt, ss, &aas, &donor, &acceptor); \
	ns_gen_prof(ns_int_t, km, aas, al, opt, neg_inf, &mem_ap, &ap); \
	go = sse_gen(set1, _suf)(opt->go); \
	ge = sse_gen(set1, _suf)(opt->ge); \
	goe= sse_gen(set1, _suf)(opt->go + opt->ge); \
	io = sse_gen(set1, _suf)(opt->io); \
	fs = sse_gen(set1, _suf)(opt->fs); \
	H0 = ns_alloc16(km, (slen + 1) * 5 + slen * 7, &mem_H); \
	H = H0 + 1, H1 = H0 + (slen + 1) + 1, H2 = H0 + (slen + 1) * 2 + 1, H3 = H0 + (slen + 1) * 3 + 1; \
	Hmax = H0 + (slen + 1) * 4 + 1; \
	D = Hmax + slen, D1 = D + slen, D2 = D1 + slen, D3 = D2 + slen; \
	A = D3 + slen, B = A + slen, C = B + slen; \
	if ((opt->flag & NS_F_CIGAR) && !is_ext) \
		tb = ns_alloc16(km, nl * slen, &mem_tb);

#define NS_GEN_INIT1(_suf) \
	for (i = 0; i < (slen + 1) * 5 + slen * 7; ++i) \
		H0[i] = sse_gen(set1, _suf)(neg_inf); \
	H3[-1] = sse_gen(insert, _suf)(H3[-1], 0, 0); \
	H2[-1] = sse_gen(insert, _suf)(H2[-1], -opt->fs, 0); \
	H1[-1] = sse_gen(insert, _suf)(H1[-1], -opt->fs, 0);

#define NS_GEN_INIT2(_suf) \
		int32_t k; \
		__m128i *tmp, I, *S = ap + nas[i] * slen, dim1, di, dip1, ai, aim1, aim2, last_h, gei; \
		gei = nas[i] == 20? fs : ge; \
		dim1 = sse_gen(set1, _suf)(donor[i-1]), di = sse_gen(set1, _suf)(donor[i]), dip1 = sse_gen(set1, _suf)(donor[i+1]); \
		ai = sse_gen(set1, _suf)(acceptor[i]), aim1 = sse_gen(set1, _suf)(acceptor[i-1]), aim2 = sse_gen(set1, _suf)(acceptor[i-2]); \
		I = last_h = sse_gen(set1, _suf)(neg_inf); \
		if (i > 2) { /* FIXME: this is close but not correct */ \
			H3[-1] = sse_gen(insert, _suf)(_mm_slli_si128(H3[slen - 1], sizeof(ns_int_t)), neg_inf, 0); \
			H2[-1] = sse_gen(insert, _suf)(_mm_slli_si128(H2[slen - 1], sizeof(ns_int_t)), neg_inf, 0); \
			H1[-1] = sse_gen(insert, _suf)(_mm_slli_si128(H1[slen - 1], sizeof(ns_int_t)), neg_inf, 0); \
		}

static inline int ns_le_epi16(__m128i a, __m128i b)
{
#if defined(__SSE2__)
	return !_mm_movemask_epi8(_mm_cmpgt_epi16(a, b));
#elif defined(__ARM_NEON)
	return (vmaxvq_u8(_mm_cmpgt_epi16(a, b)) == 0);
#endif
}

static inline int ns_le_epi32(__m128i a, __m128i b)
{
#if defined(__SSE2__)
	return !_mm_movemask_epi8(_mm_cmpgt_epi16(a, b));
#elif defined(__ARM_NEON)
	return (vmaxvq_u8(_mm_cmpgt_epi32(a, b)) == 0);
#endif
}

static inline __m128i ns_select(__m128i cond, __m128i a, __m128i b)
{
#if defined(__ARM_NEON) || defined(__SSE4_1__) // there is an ARM emulation of _mm_blendv_epi8()
	return _mm_blendv_epi8(b, a, cond);
#elif defined(__SSE2__)
	return _mm_or_si128(_mm_and_si128(a, cond), _mm_andnot_si128(cond, b));
#endif
}

static inline int ns_max_8(__m128i a)
{
#if defined(__ARM_NEON)
	return vmaxvq_s16(vreinterpretq_s16_u8(a));
#elif defined(__SSE2__)
	a = _mm_max_epi16(a, _mm_srli_si128(a, 8));
	a = _mm_max_epi16(a, _mm_srli_si128(a, 4));
	a = _mm_max_epi16(a, _mm_srli_si128(a, 2));
    return (int16_t)_mm_extract_epi16(a, 0);
#endif
}

#if defined(__SSE2__) && !defined(__SSE4_1__)
static inline __m128i _mm_max_epi32(__m128i a, __m128i b) { return ns_select(_mm_cmpgt_epi32(a, b), a, b); }
static inline __m128i _mm_insert_epi32(__m128i a, int b, const int ndx)
{
    switch (ndx & 0x3) {
    case 0: a = _mm_insert_epi16( a, b    , 0 );
            a = _mm_insert_epi16( a, b<<16, 1 ); break;
    case 1: a = _mm_insert_epi16( a, b    , 2 );
            a = _mm_insert_epi16( a, b<<16, 3 ); break;
    case 2: a = _mm_insert_epi16( a, b    , 4 );
            a = _mm_insert_epi16( a, b<<16, 5 ); break;
    case 3: a = _mm_insert_epi16( a, b    , 6 );
            a = _mm_insert_epi16( a, b<<16, 7 ); break;
    }
    return a;
}
#endif

static inline float ns_log2(float x) // NB: this doesn't work when x<2
{
	union { float f; uint32_t i; } z = { x };
	float log_2 = ((z.i >> 23) & 255) - 128;
	z.i &= ~(255 << 23);
	z.i += 127 << 23;
	log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
	return log_2;
}

void ns_global_gs16b(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, const uint8_t *ss, ns_rst_t *r)
{
	typedef int16_t ns_int_t;
	const int32_t ssize = sizeof(ns_int_t), vsize = 16 / ssize;
	const ns_int_t neg_inf = (ns_int_t)(1 << (8*ssize - 1));
	NS_GEN_VAR(epi16)
	NS_GEN_PREPARE(epi16)
	NS_GEN_INIT1(epi16)

	if (tb == 0) {
		int32_t max_sc = INT32_MIN, max_sc_log = INT32_MIN, tmp_sc, tmp_sc_log, end_sc, max_i = -1, pen_len = al * 3;
		for (i = 2; i < nl; ++i) {
			__m128i max;
			NS_GEN_INIT2(epi16)
			max = _mm_set1_epi16(neg_inf);
			for (j = 0; j < slen; ++j) {
				__m128i h, t, u, v;
				// H(i-3,j-1) + s(i,j)
				u = _mm_load_si128(H3 + j - 1);
				v = _mm_load_si128(S + j);
				h = _mm_adds_epi16(u, v);
				// I(i,j) = max{ H(i,j-1) - q, I(i,j-1) } - e
				t = _mm_subs_epi16(last_h, go);
				t = _mm_max_epi16(t, I);
				I = _mm_subs_epi16(t, ge);
				h = _mm_max_epi16(h, I);
				// D(i,j) = max{ H(i-3,j) - q, D(i-3,j) } - e
				u = _mm_load_si128(H3 + j);
				v = _mm_load_si128(D3 + j);
				t = _mm_max_epi16(_mm_subs_epi16(u, go), v);
				t = _mm_subs_epi16(t, gei);
				_mm_store_si128(D + j, t);
				h = _mm_max_epi16(h, t);
				// A(i,j) = max{ H(i-1,j)   - r - d(i-1), A(i-1,j) }
				u = _mm_subs_epi16(_mm_load_si128(H1 + j), io);
				v = _mm_load_si128(A + j);
				t = _mm_subs_epi16(u, dim1);
				t = _mm_max_epi16(t, v);
				_mm_store_si128(A + j, t);
				h = _mm_max_epi16(h, _mm_subs_epi16(t, ai));
				// B(i,j) = max{ H(i-1,j-1) - r - d(i),   B(i-1,j) }
				u = _mm_subs_epi16(_mm_load_si128(H1 + j - 1), io);
				v = _mm_load_si128(B + j);
				t = _mm_subs_epi16(u, di);
				t = _mm_max_epi16(t, v);
				_mm_store_si128(B + j, t);
				h = _mm_max_epi16(h, _mm_subs_epi16(t, aim2));
				// C(i,j) = max{ H(i-1,j-1) - r - d(i+1), C(i-1,j) }
				v = _mm_load_si128(C + j);
				t = _mm_subs_epi16(u, dip1);
				t = _mm_max_epi16(t, v);
				_mm_store_si128(C + j, t);
				h = _mm_max_epi16(h, _mm_subs_epi16(t, aim1));
				// H(i-1,j-1)-f and H(i-2,j-1)-f
				t = _mm_subs_epi16(_mm_load_si128(H1 + j), fs);
				h = _mm_max_epi16(h, t);
				t = _mm_subs_epi16(_mm_load_si128(H2 + j), fs);
				h = _mm_max_epi16(h, t);
				// H(i-1,j)-f and H(i-2,j)-f
				t = _mm_subs_epi16(_mm_load_si128(H1 + j - 1), fs);
				h = _mm_max_epi16(h, t);
				t = _mm_subs_epi16(_mm_load_si128(H2 + j - 1), fs);
				h = _mm_max_epi16(h, t);
				// save H
				max = _mm_max_epi16(max, h);
				_mm_store_si128(H + j, h);
				last_h = h;
			}
			I = _mm_max_epi16(_mm_subs_epi16(last_h, goe), _mm_subs_epi16(I, ge));
			for (k = 0; k < vsize; ++k) { // lazy-F loop
				I = _mm_insert_epi16(_mm_slli_si128(I, sizeof(ns_int_t)), neg_inf, 0);
				for (j = 0; j < slen; ++j) {
					__m128i h;
					h = _mm_load_si128(H + j);
					h = _mm_max_epi16(h, I);
					max = _mm_max_epi16(max, h);
					_mm_store_si128(H + j, h);
					h = _mm_subs_epi16(h, goe);
					I = _mm_subs_epi16(I, ge);
					if (ns_le_epi16(I, h)) break;
				}
				if (j < slen) break;
			}
			tmp_sc = ns_max_8(max);
			end_sc = *((ns_int_t*)&H[(al-1)%slen] + (al-1)/slen) + opt->end_bonus;
			tmp_sc = tmp_sc > end_sc? tmp_sc : end_sc;
			tmp_sc_log = tmp_sc - (i - pen_len < 2? 0 : (int32_t)(opt->ie_coef * ns_log2(i - pen_len) + .5f));
			if (tmp_sc_log > max_sc_log) {
				max_sc = tmp_sc, max_sc_log = tmp_sc_log, max_i = i;
				memcpy(&Hmax[-1], &H[-1], (slen + 1) * 16);
			}
			tmp = H3, H3 = H2, H2 = H1, H1 = H, H = tmp;
			tmp = D3, D3 = D2, D2 = D1, D1 = D, D = tmp;
			if (max_sc_log - tmp_sc_log > opt->xdrop) break;
		}
		if (is_ext) {
			for (j = 0; j < al; ++j) {
				int32_t sc = *((ns_int_t*)&Hmax[j%slen] + j/slen);
				if (j == al - 1) sc += opt->end_bonus;
				if (sc == max_sc) break;
			}
			assert(j < al);
			r->nt_len = max_i + 1, r->aa_len = j + 1, r->score = max_sc;
		}
	} else {
		for (i = 2; i < nl; ++i) {
			__m128i *tbi = tb + i * slen;
			NS_GEN_INIT2(epi16)
			for (j = 0; j < slen; ++j) {
				__m128i h, t, u, v, y, z;
				// H(i-3,j-1) + s(i,j)
				y = _mm_setzero_si128();
				z = _mm_setzero_si128();
				u = _mm_load_si128(H3 + j - 1);
				v = _mm_load_si128(S + j);
				h = _mm_adds_epi16(u, v);
				// I(i,j) = max{ H(i,j-1) - q, I(i,j-1) } - e
				t = _mm_subs_epi16(last_h, go);
				z = _mm_or_si128(z, _mm_and_si128(_mm_cmpgt_epi16(I, t), _mm_set1_epi16(1<<4)));
				t = _mm_max_epi16(t, I);
				I = _mm_subs_epi16(t, ge);
				y = ns_select(_mm_cmpgt_epi16(I, h), _mm_set1_epi16(1), y);
				h = _mm_max_epi16(h, I);
				// D(i,j) = max{ H(i-3,j) - q, D(i-3,j) } - e
				u = _mm_subs_epi16(_mm_load_si128(H3 + j), go);
				v = _mm_load_si128(D3 + j);
				z = _mm_or_si128(z, _mm_and_si128(_mm_cmpgt_epi16(v, u), _mm_set1_epi16(1<<5)));
				t = _mm_max_epi16(u, v);
				t = _mm_subs_epi16(t, gei);
				_mm_store_si128(D + j, t);
				y = ns_select(_mm_cmpgt_epi16(t, h), _mm_set1_epi16(2), y);
				h = _mm_max_epi16(h, t);
				// A(i,j) = max{ H(i-1,j)   - r - d(i-1), A(i-1,j) }
				u = _mm_subs_epi16(_mm_load_si128(H1 + j), io);
				v = _mm_load_si128(A + j);
				t = _mm_subs_epi16(u, dim1);
				z = _mm_or_si128(z, _mm_and_si128(_mm_cmpgt_epi16(v, t), _mm_set1_epi16(1<<6)));
				t = _mm_max_epi16(t, v);
				_mm_store_si128(A + j, t);
				t = _mm_subs_epi16(t, ai);
				y = ns_select(_mm_cmpgt_epi16(t, h), _mm_set1_epi16(3), y);
				h = _mm_max_epi16(h, t);
				// B(i,j) = max{ H(i-1,j-1) - r - d(i),   B(i-1,j) }
				u = _mm_subs_epi16(_mm_load_si128(H1 + j - 1), io);
				v = _mm_load_si128(B + j);
				t = _mm_subs_epi16(u, di);
				z = _mm_or_si128(z, _mm_and_si128(_mm_cmpgt_epi16(v, t), _mm_set1_epi16(1<<7)));
				t = _mm_max_epi16(t, v);
				_mm_store_si128(B + j, t);
				t = _mm_subs_epi16(t, aim2);
				y = ns_select(_mm_cmpgt_epi16(t, h), _mm_set1_epi16(4), y);
				h = _mm_max_epi16(h, t);
				// C(i,j) = max{ H(i-1,j-1) - r - d(i+1), C(i-1,j) }
				v = _mm_load_si128(C + j);
				t = _mm_subs_epi16(u, dip1);
				z = _mm_or_si128(z, _mm_and_si128(_mm_cmpgt_epi16(v, t), _mm_set1_epi16(1<<8)));
				t = _mm_max_epi16(t, v);
				_mm_store_si128(C + j, t);
				t = _mm_subs_epi16(t, aim1);
				y = ns_select(_mm_cmpgt_epi16(t, h), _mm_set1_epi16(5), y);
				h = _mm_max_epi16(h, t);
				// H(i-1,j-1)-f and H(i-2,j-1)-f
				t = _mm_subs_epi16(_mm_load_si128(H1 + j), fs);
				y = ns_select(_mm_cmpgt_epi16(t, h), _mm_set1_epi16(6), y);
				h = _mm_max_epi16(h, t);
				t = _mm_subs_epi16(_mm_load_si128(H2 + j), fs);
				y = ns_select(_mm_cmpgt_epi16(t, h), _mm_set1_epi16(7), y);
				h = _mm_max_epi16(h, t);
				// H(i-1,j)-f and H(i-2,j)-f
				t = _mm_subs_epi16(_mm_load_si128(H1 + j - 1), fs);
				y = ns_select(_mm_cmpgt_epi16(t, h), _mm_set1_epi16(8), y);
				h = _mm_max_epi16(h, t);
				t = _mm_subs_epi16(_mm_load_si128(H2 + j - 1), fs);
				y = ns_select(_mm_cmpgt_epi16(t, h), _mm_set1_epi16(9), y);
				h = _mm_max_epi16(h, t);
				// save H and traceback
				z = _mm_or_si128(z, y);
				_mm_store_si128(tbi + j, z);
				_mm_store_si128(H + j, h);
				last_h = h;
			}
			I = _mm_max_epi16(_mm_subs_epi16(last_h, goe), _mm_subs_epi16(I, ge));
			for (k = 0; k < vsize; ++k) { // lazy-F loop
				I = _mm_insert_epi16(_mm_slli_si128(I, sizeof(ns_int_t)), neg_inf, 0);
				for (j = 0; j < slen; ++j) {
					__m128i h, z;
					z = _mm_load_si128(tbi + j);
					h = _mm_load_si128(H + j);
					z = _mm_or_si128(z, _mm_and_si128(_mm_cmpgt_epi16(I, h), _mm_set1_epi16(1<<9)));
					h = _mm_max_epi16(h, I);
					_mm_store_si128(tbi + j, z);
					_mm_store_si128(H + j, h);
					h = _mm_subs_epi16(h, goe);
					I = _mm_subs_epi16(I, ge);
					if (ns_le_epi16(I, h)) break;
				}
				if (j < slen) break;
			}
			//printf("[%d:%d]\t", i, k); for (j = 0; j < al; ++j) printf("%d\t", *((ns_int_t*)&H[0] + j)); putchar('\n');
			tmp = H3, H3 = H2, H2 = H1, H1 = H, H = tmp;
			tmp = D3, D3 = D2, D2 = D1, D1 = D, D = tmp;
		}
	}
	if (!is_ext) r->score = *((ns_int_t*)&H1[(al-1)%slen] + (al-1)/slen);
	kfree(km, mem_H);
	kfree(km, mem_ap);
	kfree(km, nas);
	if (tb) {
		ns_backtrack(km, vsize, tb, nl, al, &r->cigar, &r->n_cigar, &r->m_cigar);
		kfree(km, mem_tb);
	}
}

void ns_global_gs16(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, ns_rst_t *r)
{
	ns_global_gs16b(km, ns, nl, as, al, opt, 0, r);
}

void ns_global_gs32b(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, const uint8_t *ss, ns_rst_t *r)
{
	typedef int32_t ns_int_t;
	const int32_t ssize = sizeof(ns_int_t), vsize = 16 / ssize;
	const ns_int_t neg_inf = -0x40000000;
	NS_GEN_VAR(epi32)
	NS_GEN_PREPARE(epi32)
	NS_GEN_INIT1(epi32)

	if (tb == 0) {
		for (i = 2; i < nl; ++i) {
			NS_GEN_INIT2(epi32)
			for (j = 0; j < slen; ++j) {
				__m128i h, t, u, v;
				// H(i-3,j-1) + s(i,j)
				u = _mm_load_si128(H3 + j - 1);
				v = _mm_load_si128(S + j);
				h = _mm_add_epi32(u, v);
				// I(i,j) = max{ H(i,j-1) - q, I(i,j-1) } - e
				t = _mm_sub_epi32(last_h, go);
				t = _mm_max_epi32(t, I);
				I = _mm_sub_epi32(t, ge);
				h = _mm_max_epi32(h, I);
				// D(i,j) = max{ H(i-3,j) - q, D(i-3,j) } - e
				u = _mm_load_si128(H3 + j);
				v = _mm_load_si128(D3 + j);
				t = _mm_max_epi32(_mm_sub_epi32(u, go), v);
				t = _mm_sub_epi32(t, gei);
				_mm_store_si128(D + j, t);
				h = _mm_max_epi32(h, t);
				// A(i,j) = max{ H(i-1,j)   - r - d(i-1), A(i-1,j) }
				u = _mm_sub_epi32(_mm_load_si128(H1 + j), io);
				v = _mm_load_si128(A + j);
				t = _mm_sub_epi32(u, dim1);
				t = _mm_max_epi32(t, v);
				_mm_store_si128(A + j, t);
				h = _mm_max_epi32(h, _mm_sub_epi32(t, ai));
				// B(i,j) = max{ H(i-1,j-1) - r - d(i),   B(i-1,j) }
				u = _mm_sub_epi32(_mm_load_si128(H1 + j - 1), io);
				v = _mm_load_si128(B + j);
				t = _mm_sub_epi32(u, di);
				t = _mm_max_epi32(t, v);
				_mm_store_si128(B + j, t);
				h = _mm_max_epi32(h, _mm_sub_epi32(t, aim2));
				// C(i,j) = max{ H(i-1,j-1) - r - d(i+1), C(i-1,j) }
				v = _mm_load_si128(C + j);
				t = _mm_sub_epi32(u, dip1);
				t = _mm_max_epi32(t, v);
				_mm_store_si128(C + j, t);
				h = _mm_max_epi32(h, _mm_sub_epi32(t, aim1));
				// H(i-1,j-1)-f and H(i-2,j-1)-f
				t = _mm_sub_epi32(_mm_load_si128(H1 + j), fs);
				h = _mm_max_epi32(h, t);
				t = _mm_sub_epi32(_mm_load_si128(H2 + j), fs);
				h = _mm_max_epi32(h, t);
				// H(i-1,j)-f and H(i-2,j)-f
				t = _mm_sub_epi32(_mm_load_si128(H1 + j - 1), fs);
				h = _mm_max_epi32(h, t);
				t = _mm_sub_epi32(_mm_load_si128(H2 + j - 1), fs);
				h = _mm_max_epi32(h, t);
				// save H
				_mm_store_si128(H + j, h);
				last_h = h;
			}
			I = _mm_max_epi32(_mm_sub_epi32(last_h, goe), _mm_sub_epi32(I, ge));
			for (k = 0; k < vsize; ++k) { // lazy-F loop
				I = _mm_insert_epi32(_mm_slli_si128(I, sizeof(ns_int_t)), neg_inf, 0);
				for (j = 0; j < slen; ++j) {
					__m128i h;
					h = _mm_load_si128(H + j);
					h = _mm_max_epi32(h, I);
					_mm_store_si128(H + j, h);
					h = _mm_sub_epi32(h, goe);
					I = _mm_sub_epi32(I, ge);
					if (ns_le_epi32(I, h)) break;
				}
				if (j < slen) break;
			}
			tmp = H3, H3 = H2, H2 = H1, H1 = H, H = tmp;
			tmp = D3, D3 = D2, D2 = D1, D1 = D, D = tmp;
		}
	} else {
		for (i = 2; i < nl; ++i) {
			__m128i *tbi = tb + i * slen;
			NS_GEN_INIT2(epi32)
			for (j = 0; j < slen; ++j) {
				__m128i h, t, u, v, y, z;
				// H(i-3,j-1) + s(i,j)
				y = _mm_setzero_si128();
				z = _mm_setzero_si128();
				u = _mm_load_si128(H3 + j - 1);
				v = _mm_load_si128(S + j);
				h = _mm_add_epi32(u, v);
				// I(i,j) = max{ H(i,j-1) - q, I(i,j-1) } - e
				t = _mm_sub_epi32(last_h, go);
				z = _mm_or_si128(z, _mm_and_si128(_mm_cmpgt_epi32(I, t), _mm_set1_epi32(1<<4)));
				t = _mm_max_epi32(t, I);
				I = _mm_sub_epi32(t, ge);
				y = ns_select(_mm_cmpgt_epi32(I, h), _mm_set1_epi32(1), y);
				h = _mm_max_epi32(h, I);
				// D(i,j) = max{ H(i-3,j) - q, D(i-3,j) } - e
				u = _mm_sub_epi32(_mm_load_si128(H3 + j), go);
				v = _mm_load_si128(D3 + j);
				z = _mm_or_si128(z, _mm_and_si128(_mm_cmpgt_epi32(v, u), _mm_set1_epi32(1<<5)));
				t = _mm_max_epi32(u, v);
				t = _mm_sub_epi32(t, gei);
				_mm_store_si128(D + j, t);
				y = ns_select(_mm_cmpgt_epi32(t, h), _mm_set1_epi32(2), y);
				h = _mm_max_epi32(h, t);
				// A(i,j) = max{ H(i-1,j)   - r - d(i-1), A(i-1,j) }
				u = _mm_sub_epi32(_mm_load_si128(H1 + j), io);
				v = _mm_load_si128(A + j);
				t = _mm_sub_epi32(u, dim1);
				z = _mm_or_si128(z, _mm_and_si128(_mm_cmpgt_epi32(v, t), _mm_set1_epi32(1<<6)));
				t = _mm_max_epi32(t, v);
				_mm_store_si128(A + j, t);
				t = _mm_sub_epi32(t, ai);
				y = ns_select(_mm_cmpgt_epi32(t, h), _mm_set1_epi32(3), y);
				h = _mm_max_epi32(h, t);
				// B(i,j) = max{ H(i-1,j-1) - r - d(i),   B(i-1,j) }
				u = _mm_sub_epi32(_mm_load_si128(H1 + j - 1), io);
				v = _mm_load_si128(B + j);
				t = _mm_sub_epi32(u, di);
				z = _mm_or_si128(z, _mm_and_si128(_mm_cmpgt_epi32(v, t), _mm_set1_epi32(1<<7)));
				t = _mm_max_epi32(t, v);
				_mm_store_si128(B + j, t);
				t = _mm_sub_epi32(t, aim2);
				y = ns_select(_mm_cmpgt_epi32(t, h), _mm_set1_epi32(4), y);
				h = _mm_max_epi32(h, t);
				// C(i,j) = max{ H(i-1,j-1) - r - d(i+1), C(i-1,j) }
				v = _mm_load_si128(C + j);
				t = _mm_sub_epi32(u, dip1);
				z = _mm_or_si128(z, _mm_and_si128(_mm_cmpgt_epi32(v, t), _mm_set1_epi32(1<<8)));
				t = _mm_max_epi32(t, v);
				_mm_store_si128(C + j, t);
				t = _mm_sub_epi32(t, aim1);
				y = ns_select(_mm_cmpgt_epi32(t, h), _mm_set1_epi32(5), y);
				h = _mm_max_epi32(h, t);
				// H(i-1,j-1)-f and H(i-2,j-1)-f
				t = _mm_sub_epi32(_mm_load_si128(H1 + j), fs);
				y = ns_select(_mm_cmpgt_epi32(t, h), _mm_set1_epi32(6), y);
				h = _mm_max_epi32(h, t);
				t = _mm_sub_epi32(_mm_load_si128(H2 + j), fs);
				y = ns_select(_mm_cmpgt_epi32(t, h), _mm_set1_epi32(7), y);
				h = _mm_max_epi32(h, t);
				// H(i-1,j)-f and H(i-2,j)-f
				t = _mm_sub_epi32(_mm_load_si128(H1 + j - 1), fs);
				y = ns_select(_mm_cmpgt_epi32(t, h), _mm_set1_epi32(8), y);
				h = _mm_max_epi32(h, t);
				t = _mm_sub_epi32(_mm_load_si128(H2 + j - 1), fs);
				y = ns_select(_mm_cmpgt_epi32(t, h), _mm_set1_epi32(9), y);
				h = _mm_max_epi32(h, t);
				// save H and traceback
				z = _mm_or_si128(z, y);
				_mm_store_si128(tbi + j, z);
				_mm_store_si128(H + j, h);
				last_h = h;
			}
			I = _mm_max_epi32(_mm_sub_epi32(last_h, goe), _mm_sub_epi32(I, ge));
			for (k = 0; k < vsize; ++k) { // lazy-F loop
				I = _mm_insert_epi32(_mm_slli_si128(I, sizeof(ns_int_t)), neg_inf, 0);
				for (j = 0; j < slen; ++j) {
					__m128i h, z;
					z = _mm_load_si128(tbi + j);
					h = _mm_load_si128(H + j);
					z = _mm_or_si128(z, _mm_and_si128(_mm_cmpgt_epi32(I, h), _mm_set1_epi32(1<<9)));
					h = _mm_max_epi32(h, I);
					_mm_store_si128(tbi + j, z);
					_mm_store_si128(H + j, h);
					h = _mm_sub_epi32(h, goe);
					I = _mm_sub_epi32(I, ge);
					if (ns_le_epi32(I, h)) break;
				}
				if (j < slen) break;
			}
			tmp = H3, H3 = H2, H2 = H1, H1 = H, H = tmp;
			tmp = D3, D3 = D2, D2 = D1, D1 = D, D = tmp;
		}
	}
	r->score = *((ns_int_t*)&H1[(al-1)%slen] + (al-1)/slen);
	kfree(km, mem_H);
	kfree(km, mem_ap);
	kfree(km, nas);
	if (tb) {
		ns_backtrack(km, vsize, tb, nl, al, &r->cigar, &r->n_cigar, &r->m_cigar);
		kfree(km, mem_tb);
	}
}

void ns_global_gs32(void *km, const char *ns, int32_t nl, const char *as, int32_t al, const ns_opt_t *opt, ns_rst_t *r)
{
	ns_global_gs32b(km, ns, nl, as, al, opt, 0, r);
}
