#include <string.h>
#include "miniprot.h"

void mp_idxopt_init(mp_idxopt_t *io)
{
	memset(io, 0, sizeof(*io));
	io->bbit = 8;
	io->min_aa_len = 15;
	io->kmer = 6;
	io->smer = 4;
}

void mp_mapopt_init(mp_mapopt_t *mo)
{
	memset(mo, 0, sizeof(*mo));
	mo->flag = 0;
	mo->mini_batch_size = 100000000;
	mo->max_occ = 50000;
	mo->max_gap = 1000;
	mo->max_intron = 200000;
	mo->bw = mo->max_intron;
	mo->min_chn_cnt = 5;
}
