Release 0.17-r279 (15 June 2025)
--------------------------------

Notable changes:

 * Bugfix: older version may output overlapping exons with `--spsc`.

(0.17: 15 June 2025, r279)



Release 0.16-r273 (22 May 2025)
-------------------------------

Notable changes:

 * Improvement: `--spsc` now works with and optimized for minisplice output.
   Note that depending on scoring in `--spsc`, miniprot may generate 1bp
   introns in extremely rare cases. This happens to 3 out of 4.6 million
   aligned introns with zebrafish minisplice scores.

This version produces alignment identical to v0.13 unless `--spsc` is used.

(0.16: 22 May 2025, r273)



Release 0.15-r270 (18 April 2025)
---------------------------------

Notable change:

 * Improvement: better alignment in tandem gene duplications (#35). This
   slightly improves overall junction accuracy on all test datasets in the
   miniprot paper. Thank the NCBI EGAPx team for providing test cases and
   explaining the relevance in gene annotation.

(0.15: 18 April 2025, r270)



Release 0.14-r265 (7 March 2025)
--------------------------------

Notable changes:

 * Bug fix: support chromosomes longer than 2Gbp (#59)

 * EXPERIMENTAL feature: read splice scores from a file specified by `--spsc`
   and consider the scores during residue alignment. The feature makes it
   possible to apply advanced splice models and to improve miniprot alignment.

 * Improvement: documented the C APIs and added an example program on using the
   APIs (#69).

This version produces alignment identical to v0.13, except for long chromosomes
or when the new feature is used.

(0.14: 7 March 2025, r265)



Release 0.13-r248 (6 March 2024)
--------------------------------

Notable changes:

 * New feature: added option -T to specify a non-standard NCBI translation
   table (#56 and #57). As this is an indexing option, the binary index format
   has to be changed accordingly. Miniprot will reject indices built with
   previous versions.

 * Improvement: properly handle reference deletions involving in-frame stop
   codons (#58). Older versions would not penalize these stop codons. This
   change also improves junction accuracy especially for distant homologs.

 * Bug fix: in the GFF3 output, CDS now includes stop codons (#55). Note the in
   GTF, CDS excludes stop codons.

 * Bug fix: suppress an extra amino acid in the --trans or --aln output (#47).
   In rare cases, this may lead to memory violation.

(0.13: 6 March 2024, r248)



Release 0.12-r237 (24 June 2023)
--------------------------------

Notable changes:

 * New feature: added option --no-cs to disable the cs tag. This tag is not as
   useful as the cs tag for nucleotide alignment because it does not encode the
   matching amino acids.

 * New feature: output the number of frameshifts and in-frame stop codons in
   the PAF output. It is non-trivial to parse in-frame stop codons.

 * Bugfix: fixed malformatted protein sequences when --gtf and --trans are both
   specified (#45).

(0.12: 24 June 2023, r237)



Release 0.11-r234 (18 April 2023)
---------------------------------

Notable changes:

 * New feature: added option --trans to output translated protein sequences. It
   is possible to extract these sequences from the --aln output but the --trans
   output is smaller and more convenient.

 * Bugfix: infinite error messages if a wrong option is in use.

 * Improvement: better error messages given nonexisting query files (#40).

(0.11: 18 April 2023, r234)



Release 0.10-r225 (3 April 2023)
--------------------------------

Notable change:

 * Bugfix: rare segmentation fault (#38 and #39). This bug affected all
   previous versions of miniprot.

(0.10: 3 April 2023, r225)



Release 0.9-r223 (9 March 2023)
-------------------------------

Notable change:

 * Bugfix: some query proteins were not outputted with option `-u`.

(0.9: 9 March 2023, r223)



Release 0.8-r220 (6 March 2023)
-------------------------------

Notable changes:

 * Improvement: slightly improved the sensitivity to distant homolog at a minor
   cost of specificity. On the human-zebrafish dataset, we gained 1.2% junction
   sensivity at the cost of 0.2% of specificity.

 * New feature: added option --aln to output residue alignment.

 * New feature: added option -I to automatically set the maximum intron size to
   sqrt(GenomeSize) * 3.6, where GenomeSize is the total length of the
   nucleotide sequences. For smaller genomes, a smaller threshold leads to
   higher accuracy. This option is not the default because the reference is not
   always a genome.

(0.8: 6 March 2023, r220)



Release 0.7-r207 (25 December 2022)
-----------------------------------

Notable changes:

 * Improvement: replaced open syncmers with modimers. This simplifies the code
   and slightly reduces the memory at comparable k-mer sampling rate. This
   changes the index format.

 * Improvement: fine tune parameters for higher sensitivity at a minor cost of
   junction accuracy: a) only index ORFs >= 30bp; b) reduced max k-mer
   occurrences from 50k to 20k; c) sample k-mers at a rate of 50%; d) reduced
   min number of k-mers from 5 to 3; e) add a bonus chaining score for anchors
   on the same reference block.

 * Improvement: adjust the max k-mer occurrence dynamically per protein.

 * Improvement: implemented 2-level chaining like minimap2 and minigraph. This
   reduces chaining time.

 * Bugfix: fixed a rare off-by-1 memory violation

 * Bugfix: fixed a memory leak

Overall, miniprot becomes faster at slightly higher peak memory usage. It is
more sensitive to distant homologs, though the junction accuracy of additional
alignment is usually lower. Also importantly, the index format of miniprot has
been changed. Miniprot will throw an error if you use miniprot with pre-built
indices generated with older versions.

(0.7: 25 December 2022, r207)



Release 0.6-r185 (12 December 2022)
-----------------------------------

Notable changes:

 * Improvement: for each protein, only output alignments close to the best
   alignment. Also added option --outs to tune the threshold.

 * New feature: output GTF with option --gtf.

(0.6: 22 December 2022, r185)



Release 0.5-r179 (17 October 2022)
----------------------------------

Notable changes:

 * Improvement: more detailed splice model considering G|GTR..YYYNYAG|. This is
   not enabled by default. Added option `-j` to change the splice model.

 * Added the miniprot preprint. Available at http://arxiv.org/abs/2210.08052

(0.5: 17 October 2022, r179)



Release 0.4-r165 (5 October 2022)
---------------------------------

This version implements a better splice model and pays a little more effort in
aligning terminal exons. It improves both sensitivity and specificity by a few
percent.

Other notable changes:

 * Breaking change: changed -C to scale the splice model

 * Bugfix: implemented option -w (#12)

 * Bugfix: reduced the indexing time for highly fragmented genomes (#12)

 * New feature: output a Rank attribute in GFF

(0.4: 5 October 2022, r165)



Release 0.3-r137 (22 September 2022)
------------------------------------

Notable changes:

 * Improvement: fine tune parameters and heuristics: higher non-GT-AG penalty,
   higher frameshift penalty, higher penalty on in-frame stop codons and a
   small penalty on long terminal introns. Miniprot is a little more sensitive
   and a little more accurate, at a minor cost of performance.

 * New feature: richer GFF output. Miniprot now reports per-exon alignment
   score, number of frameshifts, number of in-frame stop codons and
   non-canonical donor/acceptor sequences.

 * New feature: added option `--outn` to control the number of alignments per
   protein to output.

 * New feature: added option `-P` to change the ID prefix in GFF output (#6).

 * Bug fix: fixed a segmentation fault when there are no k-mers on a reference
   sequence (#4).

(0.3: 22 September 2022, r137)



Release 0.2-r116 (12 September 2022)
------------------------------------

Notable changes:

 * New feature: output GFF3 with option `--gff`. PAF alignments are embedded in
   `##PAF` lines in the output GFF3.

 * Improvement: give a bonus score -B if extension reaches the end of a protein.
   BWA-MEM and minimap2 both have this heuristic.

 * Improvement: pay more effort to the first and the last exons. This increases
   the sensitivity at the cost of performance.

 * Improvement: increased non-canonical splicing penalty -C from 6 to 11. This
   increases overall specificity.

 * Improvement: rank an alignment with the DP score disregarding introns. This
   is to reduce the effect pseudogenes. Minimap2 uses the same strategy.

 * Bug fixes: fixed incorrect CIGAR in corner cases and patched a minor memory
   leak.

This version works on several real datasets to decent accuracy without crashing
or memory leaks. It is ready for more users.

(0.2: 12 September 2022, r116)



Release 0.1-r97 (9 September 2022)
----------------------------------

This is the first public release of miniprot, a mapper for aligning proteins to
a large genome. This release has a few issues and is generally not recommended
for production uses.

(0.1: 9 September 2022, r97)
