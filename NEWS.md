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
