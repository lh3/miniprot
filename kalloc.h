#ifndef _KALLOC_H_
#define _KALLOC_H_

#include <stddef.h> /* for size_t */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	size_t capacity, available, n_blocks, n_cores, largest;
} km_stat_t;

void *kmalloc(void *km, size_t size);
void *krealloc(void *km, void *ptr, size_t size);
void *krelocate(void *km, void *ap, size_t n_bytes);
void *kcalloc(void *km, size_t count, size_t size);
void kfree(void *km, void *ptr);

void *km_init(void);
void *km_init2(void *km_par, size_t min_core_size);
void km_destroy(void *km);
void km_stat(const void *_km, km_stat_t *s);
void km_stat_print(const void *km);

#ifdef __cplusplus
}
#endif

#define Kmalloc(km, type, cnt)       ((type*)kmalloc((km), (cnt) * sizeof(type)))
#define Kcalloc(km, type, cnt)       ((type*)kcalloc((km), (cnt), sizeof(type)))
#define Krealloc(km, type, ptr, cnt) ((type*)krealloc((km), (ptr), (cnt) * sizeof(type)))

#define Kgrow(km, type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = Krealloc(km, type, ptr, (__m)); \
		} \
	} while (0)

#define KMALLOC(km, ptr, len) ((ptr) = (__typeof__(ptr))kmalloc((km), (len) * sizeof(*(ptr))))
#define KCALLOC(km, ptr, len) ((ptr) = (__typeof__(ptr))kcalloc((km), (len), sizeof(*(ptr))))

#endif
