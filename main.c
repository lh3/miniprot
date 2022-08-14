#include <stdio.h>
#include <stdlib.h>
#include "miniprot.h"
#include "ketopt.h"

int main(int argc, char *argv[])
{
	int32_t c;
	ketopt_t o = KETOPT_INIT;
	mp_idxopt_t io;
	mp_idx_t *mi;

	mp_start();
	mp_idxopt_init(&io);
	while ((c = ketopt(&o, argc, argv, 1, "k:s:b:", 0)) >= 0) {
		if (c == 'k') io.kmer = atoi(o.arg);
		else if (c == 's') io.smer = atoi(o.arg);
		else if (c == 'b') io.bbit = atoi(o.arg);
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: miniprot [options] <ref.fa>\n");
		return 1;
	}
	mi = mp_index(&io, argv[o.ind]);

	int32_t cid = 0;
	int64_t len;
	uint8_t *seq;
	mp64_v a = {0,0,0};
	seq = (uint8_t*)malloc(mi->nt->ctg[cid].len);
	len = mp_ntseq_get(mi->nt, cid, 0, -1, 0, seq);
	mp_idx_proc_seq(0, len, seq, io.min_aa_len, io.kmer, io.smer, io.bbit, 0, &a);
	return 0;
}
