#ifndef MPPRIV_H
#define MPPRIV_H

#include <stdio.h>
#include "miniprot.h"

#ifndef kroundup64
#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct mp_bseq_file_s;
typedef struct mp_bseq_file_s mp_bseq_file_t;

typedef struct {
	int32_t l_seq;
	char *name, *seq, *comment;
} mp_bseq1_t;

// from bseq.c
mp_bseq_file_t *mp_bseq_open(const char *fn);
void mp_bseq_close(mp_bseq_file_t *fp);
mp_bseq1_t *mp_bseq_read(mp_bseq_file_t *fp, int64_t chunk_size, int with_comment, int *n_);

// from misc.c
char *mp_strdup(const char *src);
void radix_sort_mp64(uint64_t *st, uint64_t *en);

// from sys.c
double mp_realtime(void);
double mp_cputime(void);
long mp_peakrss(void);

// from ntseq.c
mp_ntdb_t *mp_ntseq_read(const char *fn);
void mp_ntseq_destroy(mp_ntdb_t *db);
void mp_ntseq_dump(FILE *fp, const mp_ntdb_t *nt);
mp_ntdb_t *mp_ntseq_restore(FILE *fp);
int64_t mp_ntseq_get(const mp_ntdb_t *db, int32_t cid, int64_t st, int64_t en, int32_t rev, uint8_t *seq);

static inline uint32_t mp_hash32_mask(uint32_t key, uint32_t mask)
{
	key  = (key + ~(key << 15)) & mask;
	key ^=  (key >> 10);
	key  =  (key + (key <<  3)) & mask;
	key ^=  (key >> 6);
	key  = (key + ~(key << 11)) & mask;
	key ^=  (key >> 16);
	return key;
}

#ifdef __cplusplus
}
#endif

#endif
