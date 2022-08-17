#include <stdio.h>
#include <stdlib.h>
#include "mppriv.h"
#include "ketopt.h"

static ko_longopt_t long_options[] = {
	{ "dbg-qname",       ko_no_argument,       501 },
	{ 0, 0, 0 }
};

static inline int64_t mp_parse_num2(const char *str, char **q)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9, ++p;
	else if (*p == 'M' || *p == 'm') x *= 1e6, ++p;
	else if (*p == 'K' || *p == 'k') x *= 1e3, ++p;
	if (q) *q = p;
	return (int64_t)(x + .499);
}

static inline int64_t mp_parse_num(const char *str)
{
	return mp_parse_num2(str, 0);
}

static void print_usage(FILE *fp, const mp_idxopt_t *io, const mp_mapopt_t *mo, int n_threads)
{
	fprintf(fp, "Usage: miniprot [options] <ref.fa> <query.faa> [...]\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  Indexing:\n");
	fprintf(fp, "    -k INT       k-mer size [%d]\n", io->kmer);
	fprintf(fp, "    -s INT       submer size (density: 1/(2*(k-s)+1)) [%d]\n", io->smer);
	fprintf(fp, "    -b INT       bits per block [%d]\n", io->bbit);
	fprintf(fp, "    -d FILE      save index to FILE []\n");
	fprintf(fp, "  Mapping:\n");
	fprintf(fp, "    -c NUM       max k-mer occurrence [%d]\n", mo->max_occ);
	fprintf(fp, "    -n NUM       minimum number of syncmers in a chain [%d]\n", mo->min_chn_cnt);
	fprintf(fp, "  Input/output:\n");
	fprintf(fp, "    -t INT       number of threads [%d]\n", n_threads);
	fprintf(fp, "    -K NUM       query batch size [%ld]\n", (long)mo->mini_batch_size);
}

int main(int argc, char *argv[])
{
	int32_t c, i, n_threads = 4;
	ketopt_t o = KETOPT_INIT;
	mp_mapopt_t mo;
	mp_idxopt_t io;
	mp_idx_t *mi;
	char *fn_idx = 0;

	mp_start();
	mp_mapopt_init(&mo);
	mp_idxopt_init(&io);
	while ((c = ketopt(&o, argc, argv, 1, "k:s:b:t:d:c:n:K:", long_options)) >= 0) {
		if (c == 'k') io.kmer = atoi(o.arg);
		else if (c == 's') io.smer = atoi(o.arg);
		else if (c == 'b') io.bbit = atoi(o.arg);
		else if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'd') fn_idx = o.arg;
		else if (c == 'c') mo.max_occ = mp_parse_num(o.arg);
		else if (c == 'n') mo.min_chn_cnt = mp_parse_num(o.arg);
		else if (c == 'K') mo.mini_batch_size = mp_parse_num(o.arg);
		else if (c == 501) mp_dbg_flag |= MP_DBG_QNAME;
	}
	if (argc - o.ind == 0 || (argc - o.ind == 1 && fn_idx == 0)) {
		print_usage(stderr, &io, &mo, n_threads);
		return 1;
	}

	mi = mp_idx_load(argv[o.ind], &io, n_threads);
	if (mi == 0) {
		if (mp_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m failed to open/build the index\033[0m\n");
		return 1;
	}
	if (mp_verbose >= 3) mp_idx_print_stat(mi, mo.max_occ);
	if (fn_idx != 0) mp_idx_dump(fn_idx, mi);
	for (i = o.ind + 1; i < argc; ++i)
		mp_map_file(mi, argv[i], &mo, n_threads);
	mp_idx_destroy(mi);

	if (mp_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MP_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, mp_realtime(), mp_cputime(), mp_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
