#include <stdio.h>
#include <stdlib.h>
#include "miniprot.h"
#include "ketopt.h"

static int usage(FILE *fp, const mp_idxopt_t *io, int n_threads)
{
	fprintf(fp, "Usage: miniprot [options] <ref.fa>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  Indexing:\n");
	fprintf(fp, "    -k INT       k-mer size [%d]\n", io->kmer);
	fprintf(fp, "    -s INT       submer size [%d]\n", io->smer);
	fprintf(fp, "    -b INT       bits per block [%d]\n", io->bbit);
	fprintf(fp, "    -d FILE      save index to FILE []\n");
	fprintf(fp, "  Input/output:\n");
	fprintf(fp, "    -t INT       number of threads [%d]\n", n_threads);
	return 1;
}

int main(int argc, char *argv[])
{
	int32_t c, n_threads = 4;
	ketopt_t o = KETOPT_INIT;
	mp_idxopt_t io;
	mp_idx_t *mi;
	char *fn_idx = 0;

	mp_start();
	mp_idxopt_init(&io);
	while ((c = ketopt(&o, argc, argv, 1, "k:s:b:t:d:", 0)) >= 0) {
		if (c == 'k') io.kmer = atoi(o.arg);
		else if (c == 's') io.smer = atoi(o.arg);
		else if (c == 'b') io.bbit = atoi(o.arg);
		else if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'd') fn_idx = o.arg;
	}
	if (argc - o.ind < 1)
		return usage(stderr, &io, n_threads);
	mi = mp_idx_read(argv[o.ind], &io, n_threads);
	if (fn_idx != 0) mp_idx_dump(fn_idx, mi);

	mp_idx_destroy(mi);
	return 0;
}
