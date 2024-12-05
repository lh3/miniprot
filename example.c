// To compile:
//   gcc -g -O2 example.c libminiprot.a -lz

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "miniprot.h"
#include "nasw.h" // needed for CIGAR
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	mp_idxopt_t iopt;
	mp_mapopt_t mopt;
	int n_threads = 3;

	mp_start(); // this set the translation table to "1"
	//ns_make_tables(1); // set the translation table if not "1"
	mp_verbose = 2; // disable message output to stderr
	mp_idxopt_init(&iopt);
	mp_mapopt_init(&mopt);

	if (argc < 3) {
		fprintf(stderr, "Usage: miniprot-lite <target.fa> <query.fa>\n");
		return 1;
	}

	// open query file for reading; you may use your favorite FASTA/Q parser
	gzFile f = gzopen(argv[2], "r");
	assert(f);
	kseq_t *ks = kseq_init(f);

	// open index reader
	mp_idx_t *mi = mp_idx_load(argv[1], &iopt, n_threads);
	if (mi != 0) {
		mp_tbuf_t *tbuf = mp_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
		while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
			mp_reg1_t *reg;
			int j, i, n_reg;
			reg = mp_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query
			for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
				mp_reg1_t *r = &reg[j];
				assert(r->p);
				const mp_ctg_t *ctg = &mi->nt->ctg[r->vid>>1];
				int64_t rs = r->vid&1? ctg->len - r->ve : r->vs;
				int64_t re = r->vid&1? ctg->len - r->vs : r->ve;
				printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, (int)ks->seq.l, r->qs, r->qe, "+-"[r->vid&1]);
				printf("%s\t%ld\t%ld\t%ld\t%d\t%d\t0\tcg:Z:", ctg->name, (long)ctg->len, (long)rs, (long)re, r->p->n_iden * 3, r->p->blen);
				for (i = 0; i < r->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
					printf("%d%c", r->p->cigar[i]>>4, NS_CIGAR_STR[r->p->cigar[i]&0xf]);
				putchar('\n');
				free(r->p);
			}
			free(reg);
		}
		mp_tbuf_destroy(tbuf);
		mp_idx_destroy(mi);
	}
	kseq_destroy(ks); // close the query file
	gzclose(f);
	return 0;
}
