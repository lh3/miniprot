#include <stdio.h>
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
	int64_t i, l;
	uint8_t seq[100];
	l = mp_ntseq_get(db, 1, 1000, 1010, 0, seq);
	for (i = 0; i < l; ++i) putchar("ACGTN"[seq[i]]); putchar('\n');
	return 0;
}
