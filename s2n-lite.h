#ifndef S2N_LITE_H
#define S2N_LITE_H

#include <arm_neon.h>

typedef uint8x16_t __m128i;

static inline __m128i _mm_load_si128(const __m128i *ptr) { return vld1q_u8((const uint8_t*)ptr); }
static inline void _mm_store_si128(__m128i *ptr, __m128i a) { vst1q_u8((uint8_t*)ptr, a); }
static inline __m128i _mm_setzero_si128(void) { return vdupq_n_u8(0); }
static inline __m128i _mm_or_si128(__m128i a, __m128i b) { return vreinterpretq_u8_s32(vorrq_s32(vreinterpretq_s32_u8(a), vreinterpretq_s32_u8(b))); }
static inline __m128i _mm_and_si128(__m128i a, __m128i b) { return vreinterpretq_u8_s32(vandq_s32(vreinterpretq_s32_u8(a), vreinterpretq_s32_u8(b))); }

#define _mm_slli_si128(a, imm8) vextq_u8(_mm_setzero_si128(), (a), 16 - (imm8))

static inline __m128i _mm_blendv_epi8(__m128i a, __m128i b, __m128i mask) { return vbslq_u8(vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(mask), 7)), b, a); }

static inline __m128i _mm_set1_epi16(int a) { return vreinterpretq_u8_s16(vdupq_n_s16(a)); }
static inline __m128i _mm_cmpgt_epi16(__m128i a, __m128i b) { return vreinterpretq_u8_u16(vcgtq_s16(vreinterpretq_s16_u8(a), vreinterpretq_s16_u8(b))); }
static inline __m128i _mm_max_epi16(__m128i a, __m128i b) { return vreinterpretq_u8_s16(vmaxq_s16(vreinterpretq_s16_u8(a), vreinterpretq_s16_u8(b))); }
static inline __m128i _mm_adds_epi16(__m128i a, __m128i b) { return vreinterpretq_u8_s16(vqaddq_s16(vreinterpretq_s16_u8(a), vreinterpretq_s16_u8(b))); }
static inline __m128i _mm_subs_epi16(__m128i a, __m128i b) { return vreinterpretq_u8_s16(vqsubq_s16(vreinterpretq_s16_u8(a), vreinterpretq_s16_u8(b))); }

#define _mm_insert_epi16(a, b, imm8) vreinterpretq_u8_s16(vsetq_lane_s16((b), vreinterpretq_s16_u8(a), (imm8)))

static inline __m128i _mm_set1_epi32(int a) { return vreinterpretq_u8_s32(vdupq_n_s32(a)); }
static inline __m128i _mm_cmpgt_epi32(__m128i a, __m128i b) { return vreinterpretq_u8_u32(vcgtq_s32(vreinterpretq_s32_u8(a), vreinterpretq_s32_u8(b))); }
static inline __m128i _mm_max_epi32(__m128i a, __m128i b) { return vreinterpretq_u8_s32(vmaxq_s32(vreinterpretq_s32_u8(a), vreinterpretq_s32_u8(b))); }
static inline __m128i _mm_add_epi32(__m128i a, __m128i b) { return vreinterpretq_u8_s32(vaddq_s32(vreinterpretq_s32_u8(a), vreinterpretq_s32_u8(b))); }
static inline __m128i _mm_sub_epi32(__m128i a, __m128i b) { return vreinterpretq_u8_s32(vsubq_s32(vreinterpretq_s32_u8(a), vreinterpretq_s32_u8(b))); }

#define _mm_insert_epi32(a, b, imm8) vreinterpretq_u8_s32(vsetq_lane_s32((b), vreinterpretq_s32_u8(a), (imm8)))

#endif
