#include <stdio.h>
#include <stdlib.h>
#include "mppriv.h"
#include "ketopt.h"

static ko_longopt_t long_options[] = {
	{ "gff",             ko_no_argument,       301 },
	{ "xdrop",           ko_required_argument, 302 },
	{ "outn",            ko_required_argument, 303 },
	{ "gff-only",        ko_no_argument,       304 },
	{ "gff-delim",       ko_required_argument, 305 },
	{ "J2",              ko_required_argument, 306 },
	{ "gtf",             ko_no_argument,       307 },
	{ "outs",            ko_required_argument, 308 },
	{ "max-skip",        ko_required_argument, 309 },
	{ "no-pre-chain",    ko_no_argument,       310 },
	{ "aln",             ko_no_argument,       311 },
	{ "max-intron-out",  ko_required_argument, 312 },
	{ "outc",            ko_required_argument, 313 },
	{ "ie-coef",         ko_required_argument, 314 },
	{ "trans",           ko_no_argument,       315 },
	{ "no-cs",           ko_no_argument,       316 },
	{ "version",         ko_no_argument,       401 },
	{ "no-kalloc",       ko_no_argument,       501 },
	{ "dbg-qname",       ko_no_argument,       502 },
	{ "dbg-no-refine",   ko_no_argument,       503 },
	{ "dbg-aflt",        ko_no_argument,       504 },
	{ "dbg-anchor",      ko_no_argument,       505 },
	{ "dbg-chain",       ko_no_argument,       506 },
    { "gc",        ko_required_argument, 601 },
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
	fprintf(fp, "    -M INT       modimisers bit (sample rate = 1/2**M) [%d]\n", io->mod_bit);
	fprintf(fp, "    -L INT       min ORF length to index [%d]\n", io->min_aa_len);
	fprintf(fp, "    -b INT       bits per block [%d]\n", io->bbit);
	fprintf(fp, "    -d FILE      save index to FILE []\n");
	fprintf(fp, "  Mapping:\n");
	fprintf(fp, "    -S           no splicing (applying -G1k -J1k -e1k)\n");
	fprintf(fp, "    -c NUM       max k-mer occurrence [%d]\n", mo->max_occ);
	fprintf(fp, "    -G NUM       max intron size; override -I [200k]\n");
	fprintf(fp, "    -I           set max intron size to 3.6*sqrt(refLen)\n");
	fprintf(fp, "    -w FLOAT     weight of log gap penalty [%g]\n", mo->chn_coef_log);
	fprintf(fp, "    -n NUM       minimum number of syncmers in a chain [%d]\n", mo->min_chn_cnt);
	fprintf(fp, "    -m NUM       min chaining score [%d]\n", mo->min_chn_sc);
	fprintf(fp, "    -l INT       k-mer size for the second round of chaining [%d]\n", mo->kmer2);
	fprintf(fp, "    -e NUM       max extension for 2nd round of chaining and alignment [%d]\n", mo->max_ext);
	fprintf(fp, "    -p FLOAT     min secondary-to-primary score ratio [%g]\n", mo->pri_ratio);
	fprintf(fp, "    -N NUM       consider at most INT secondary alignments [%d]\n", mo->best_n);
	fprintf(fp, "  Alignment:\n");
//	fprintf(fp, "    -A           no alignment\n");
	fprintf(fp, "    -O INT       gap open penalty [%d]\n", mo->go);
	fprintf(fp, "    -E INT       gap extension (a k-long gap costs O+k*E) [%d]\n", mo->ge);
	fprintf(fp, "    -J INT       intron open penalty [%d]\n", mo->io);
	fprintf(fp, "    -F INT       penalty for frameshifts or in-frame stop codons [%d]\n", mo->fs);
	fprintf(fp, "    -C FLOAT     weight of splice penalty; 0 to ignore splice signals [%g]\n", mo->sp_scale);
	fprintf(fp, "    -B INT       bonus score for alignment reaching query ends [%d]\n", mo->end_bonus);
	fprintf(fp, "    -j INT       splice model: 2=mammal, 1=general, 0=none (see manual) [%d]\n", mo->sp_model);
	fprintf(fp, "    --gc=INT     genetic code: 0,1=standard, 4=mold [%d]\n", mo->gen_code);
	fprintf(fp, "  Input/output:\n");
	fprintf(fp, "    -t INT       number of threads [%d]\n", n_threads);
	fprintf(fp, "    --gff        output in the GFF3 format\n");
	fprintf(fp, "    --gtf        basic GTF output without detailed alignment\n");
	fprintf(fp, "    --aln        output residue alignment\n");
	fprintf(fp, "    --trans      output translated protein sequences (skipping frameshift)\n");
	fprintf(fp, "    -P STR       prefix for IDs in GFF3 [%s]\n", mo->gff_prefix);
	fprintf(fp, "    -u           print unmapped query proteins in PAF\n");
	fprintf(fp, "    --outn=NUM   output up to min{NUM,-N} alignments per query [%d]\n", mo->out_n);
	fprintf(fp, "    --outs=FLOAT output if score at least FLOAT*bestScore [%g]\n", mo->out_sim);
	fprintf(fp, "    --outc=FLOAT output if at least FLOAT fraction of query is aligned [%g]\n", mo->out_cov);
	fprintf(fp, "    -K NUM       query batch size [2M]\n");
}

int main(int argc, char *argv[])
{
	int32_t c, i, set_I = 0, set_G = 0, n_threads = 4;
	ketopt_t o = KETOPT_INIT;
	mp_mapopt_t mo;
	mp_idxopt_t io;
	mp_idx_t *mi;
	char *fn_idx = 0;

	mp_mapopt_init(&mo);
	mp_idxopt_init(&io);
	while ((c = ketopt(&o, argc, argv, 1, "k:M:L:s:l:b:t:d:c:n:m:K:p:N:SAO:E:J:C:F:G:e:uB:P:w:j:g:I", long_options)) >= 0) {
		if (c == 'k') io.kmer = atoi(o.arg);
		else if (c == 'M') io.mod_bit = atoi(o.arg);
		else if (c == 'L') io.min_aa_len = atoi(o.arg);
		else if (c == 'b') io.bbit = atoi(o.arg);
		else if (c == 'd') fn_idx = o.arg;
		else if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'l') mo.kmer2 = atoi(o.arg);
		else if (c == 'c') mo.max_occ = mp_parse_num(o.arg);
		else if (c == 'G') mo.bw = mo.max_intron = mp_parse_num(o.arg), set_G = 1;
		else if (c == 'I') set_I = 1;
		else if (c == 'n') mo.min_chn_cnt = mp_parse_num(o.arg);
		else if (c == 'm') mo.min_chn_sc = mp_parse_num(o.arg);
		else if (c == 'K') mo.mini_batch_size = mp_parse_num(o.arg);
		else if (c == 'p') mo.pri_ratio = atof(o.arg);
		else if (c == 'N') mo.best_n = mp_parse_num(o.arg);
		else if (c == 'S') mo.flag |= MP_F_NO_SPLICE, mo.bw = mo.max_intron = mo.max_ext = 1000, mo.io = mo.io_end = 10000, set_G = 1;
		else if (c == 'A') mo.flag |= MP_F_NO_ALIGN;
		else if (c == 'O') mo.go = atoi(o.arg);
		else if (c == 'E') mo.ge = atoi(o.arg);
		else if (c == 'J') mo.io = atoi(o.arg);
		else if (c == 'C') mo.sp_scale = atof(o.arg);
		else if (c == 'F') mp_mapopt_set_fs(&mo, atoi(o.arg));
		else if (c == 'B') mo.end_bonus = atoi(o.arg);
		else if (c == 'e') mo.max_ext = mp_parse_num(o.arg);
		else if (c == 'P') mo.gff_prefix = o.arg;
		else if (c == 'u') mo.flag |= MP_F_SHOW_UNMAP;
		else if (c == 'w') mo.chn_coef_log = atof(o.arg);
		else if (c == 'j') mo.sp_model = atoi(o.arg);
		else if (c == 'g') mo.max_gap = mp_parse_num(o.arg);
		else if (c == 301) mo.flag |= MP_F_GFF; // --gff
		else if (c == 302) mo.xdrop = atoi(o.arg); // --xdrop
		else if (c == 303) mo.out_n = mp_parse_num(o.arg); // --outn
		else if (c == 308) mo.out_sim = atof(o.arg); // --outs
		else if (c == 304) mo.flag |= MP_F_GFF | MP_F_NO_PAF; // --gff-only
		else if (c == 305) mo.gff_delim = o.arg[0]; // --gff-delim
		else if (c == 306) mo.io_end = atoi(o.arg); // --J2
		else if (c == 307) mo.flag |= MP_F_GTF; // --gtf
		else if (c == 309) mo.max_chn_max_skip = mp_parse_num(o.arg); // --max-skip
		else if (c == 310) mo.flag |= MP_F_NO_PRE_CHAIN; // --no-pre-chain
		else if (c == 311) mo.flag |= MP_F_SHOW_RESIDUE; // --aln
		else if (c == 312) mo.max_intron_flank = (mp_parse_num(o.arg) + 1) / 2; // --max-intron-out
		else if (c == 313) mo.out_cov = atof(o.arg); // --outc
		else if (c == 314) mo.ie_coef = atof(o.arg); // --ie-coef
		else if (c == 315) mo.flag |= MP_F_SHOW_TRANS; // --trans
		else if (c == 316) mo.flag |= MP_F_NO_CS; // --no-cs
		else if (c == 501) mp_dbg_flag |= MP_DBG_NO_KALLOC; // --no-kalloc
		else if (c == 502) mp_dbg_flag |= MP_DBG_QNAME; // --dbg-qname
		else if (c == 503) mp_dbg_flag |= MP_DBG_NO_REFINE; // --dbg-no-refine
		else if (c == 504) mp_dbg_flag |= MP_DBG_MORE_DP; // --dbg-aflt
		else if (c == 505) mp_dbg_flag |= MP_DBG_ANCHOR; // --dbg-anchor
		else if (c == 506) mp_dbg_flag |= MP_DBG_CHAIN; // --dbg-chain
		else if (c == 601) mo.gen_code = atoi(o.arg); // --gc
		else if (c == 's') {
			fprintf(stderr, "Option '-s' is deprecated.\n");
		} else if (c == 401) {
			printf("%s\n", MP_VERSION);
			return 0;
		} else {
			fprintf(stderr, "[WARNING]\033[1;31m unrecognized option: %s\033[0m\n", argv[o.i-1]);
		}
	}
    int codon_type = 0;
    if(mo.gen_code>0) {
        codon_type = mo.gen_code-1;
    }
	mp_start(codon_type);
	if (mp_mapopt_check(&mo) < 0) return 1;
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
	if (set_I && !set_G) mp_mapopt_set_max_intron(&mo, mi->nt->l_seq);
	if (mp_verbose >= 3) mp_idx_print_stat(mi, mo.max_occ);
	if (fn_idx != 0) mp_idx_dump(fn_idx, mi);
	for (i = o.ind + 1; i < argc; ++i) {
		int32_t res = mp_map_file(mi, argv[i], &mo, n_threads);
		if (res != 0) {
			fprintf(stderr, "[M::%s] ERROR during mapping %s (check files exists and are amino acid fastas)\n", __func__, argv[i]);
			return 1;
		}
	}
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
