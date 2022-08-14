#include <stdio.h>
#include <stdlib.h>
#include "miniprot.h"
#include "ketopt.h"

static int usage(FILE *fp, const mp_idxopt_t *io)
{
	fprintf(fp, "Usage: miniprot [options] <ref.fa>\n");
	return 1;
}

int main(int argc, char *argv[])
{
	int32_t c, n_threads = 4;
	ketopt_t o = KETOPT_INIT;
	mp_idxopt_t io;
	mp_idx_t *mi;

	mp_start();
	mp_idxopt_init(&io);
	while ((c = ketopt(&o, argc, argv, 1, "k:s:b:t:", 0)) >= 0) {
		if (c == 'k') io.kmer = atoi(o.arg);
		else if (c == 's') io.smer = atoi(o.arg);
		else if (c == 'b') io.bbit = atoi(o.arg);
		else if (c == 't') n_threads = atoi(o.arg);
	}
	if (argc - o.ind < 1)
		return usage(stderr, &io);
	mi = mp_idx_build(&io, argv[o.ind], n_threads);

	mp_idx_destroy(mi);
	return 0;
}
