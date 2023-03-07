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
