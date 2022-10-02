#include <string.h>
#include <assert.h>
#include "miniprot.h"
#include "nasw.h"

void mp_idxopt_init(mp_idxopt_t *io)
{
	memset(io, 0, sizeof(*io));
	io->bbit = 8;
	io->min_aa_len = 15;
	io->kmer = 6;
	io->smer = 4;
}

void mp_mapopt_set_fs(mp_mapopt_t *mo, int32_t fs)
{
	assert(fs >= INT8_MIN && fs <= INT8_MAX);
	mo->fs = fs;
	ns_set_stop_sc(mo->asize, mo->mat, mo->fs);
}

void mp_mapopt_init(mp_mapopt_t *mo)
{
	memset(mo, 0, sizeof(*mo));
	mo->flag = 0;
	mo->mini_batch_size = 2000000;
	mo->max_occ = 50000;
	mo->max_gap = 1000;
	mo->max_intron = 200000;
	mo->bw = mo->max_intron;
	mo->min_chn_cnt = 5;
	mo->max_chn_max_skip = 25;
	mo->max_chn_iter = 1000000;
	mo->min_chn_sc = 0;
	mo->chn_coef_log = 0.75f;
	mo->max_ext = 10000;
	mo->max_ava = 1000;
	mo->mask_level = 0.5f;
	mo->mask_len = INT32_MAX;
	mo->pri_ratio = 0.7f;
	mo->best_n = 30;
	mo->out_n = 1000;
	mo->kmer2 = 5;

	mo->go = 11, mo->ge = 1;
	mo->io = 29;
	mo->nc = 11;
	mo->fs = 17;
	mo->io_end = 19;
	mo->end_bonus = 5;
	mo->xdrop = 100;
	mo->asize = 22;
	memcpy(mo->mat, ns_mat_blosum62, 484);
	ns_set_stop_sc(mo->asize, mo->mat, mo->fs);

	mo->gff_delim = -1;
	mo->gff_prefix = "MP";
}
