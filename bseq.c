#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "kvec-km.h"
#include "mppriv.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

struct mp_bseq_file_s {
	gzFile fp;
	kseq_t *ks;
	mp_bseq1_t s;
};

mp_bseq_file_t *mp_bseq_open(const char *fn)
{
	mp_bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (f == 0) return 0;
	fp = (mp_bseq_file_t*)calloc(1, sizeof(mp_bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

void mp_bseq_close(mp_bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

static inline char *kstrdup(const kstring_t *s)
{
	char *t;
	t = (char*)malloc(s->l + 1);
	memcpy(t, s->s, s->l + 1);
	return t;
}

static inline void kseq2bseq(kseq_t *ks, mp_bseq1_t *s, int with_comment)
{
	if (ks->name.l == 0)
		fprintf(stderr, "[WARNING]\033[1;31m empty sequence name in the input.\033[0m\n");
	s->name = kstrdup(&ks->name);
	s->seq = kstrdup(&ks->seq);
	s->comment = with_comment && ks->comment.l? kstrdup(&ks->comment) : 0;
	s->l_seq = ks->seq.l;
}

mp_bseq1_t *mp_bseq_read(mp_bseq_file_t *fp, int64_t chunk_size, int with_comment, int *n_)
{
	int64_t size = 0;
	int ret;
	kvec_t(mp_bseq1_t) a = {0,0,0};
	kseq_t *ks = fp->ks;
	*n_ = 0;
	if (fp->s.seq) {
		kv_resize(mp_bseq1_t, 0, a, 256);
		kv_push(mp_bseq1_t, 0, a, fp->s);
		size = fp->s.l_seq;
		memset(&fp->s, 0, sizeof(mp_bseq1_t));
	}
	while ((ret = kseq_read(ks)) >= 0) {
		mp_bseq1_t *s;
		assert(ks->seq.l <= INT32_MAX);
		if (a.m == 0) kv_resize(mp_bseq1_t, 0, a, 256);
		kv_pushp(mp_bseq1_t, 0, a, &s);
		kseq2bseq(ks, s, with_comment);
		size += s->l_seq;
		if (size >= chunk_size)
			break;
	}
	if (ret < -1) {
		if (a.n) fprintf(stderr, "[WARNING]\033[1;31m failed to parse the FASTA/FASTQ record next to '%s'. Continue anyway.\033[0m\n", a.a[a.n-1].name);
		else fprintf(stderr, "[WARNING]\033[1;31m failed to parse the first FASTA/FASTQ record. Continue anyway.\033[0m\n");
	}
	*n_ = a.n;
	return a.a;
}
