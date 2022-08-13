#include <stdio.h>
#include <stdlib.h>
#include "miniprot.h"
#include "ketopt.h"

int main(int argc, char *argv[])
{
	int32_t c;
	ketopt_t o = KETOPT_INIT;
	mp_ntdb_t *db;

	while ((c = ketopt(&o, argc, argv, 1, "", 0)) >= 0) {
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: miniprot [options] <ref.fa>\n");
		return 1;
	}
	mp_make_tables(0);
	db = mp_ntseq_read(argv[o.ind]);

	int32_t cid = 0;
	int64_t len;
	uint8_t *seq;
	seq = (uint8_t*)malloc(db->ctg[cid].len);
	len = mp_ntseq_get(db, cid, 0, -1, 0, seq);
	mp_ntseq_get_orf(len, seq);
	return 0;
}
