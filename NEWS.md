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
