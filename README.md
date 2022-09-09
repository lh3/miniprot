**WARNING: miniprot is WIP.**

## Getting Started
```sh
# download and compile
git clone https://github.com/lh3/miniprot
cd miniprot && make

# alignment
./miniprot -ut16 test/DPP3-hs.gen.fa.gz test/DPP3-mm.pep.fa.gz > aln.paf
paftools.js paf2gff -a aln.paf > aln.gff  # requiring latest paftools.js from github

# output format
man ./miniprot.1
```

## Introduction

Miniprot aligns a protein sequence against a genome with affine gap penalty,
splicing and frameshift. It is primarily intended for annotating protein-coding
genes in a new species using known genes from other species. Miniprot is
similar to GeneWise and Exonerate in functionality but it can map proteins to
whole genomes and is much faster at the residue alignment step.

Miniprot is not optimized for mapping distant homologs because distant homologs
are less informative to gene annotations. Nonetheless, you can tune seeding
parameters to achieve higher sensitivity at the cost of performance.

## Users' Guide

### Installation

Miniprot requires SSE2 or NEON instructions and only works on x86\_64 or ARM
CPUs. It depends on zlib for parsing gzip'd input files. To compile miniprot,
type `make` in the source code directory. This will produce a standalone
executable `miniprot`. This executable is all you need to invoke miniprot.

### Usage

To run miniprot, use
```sh
miniprot -t8 ref-file protein.faa > output.paf
```
where `ref-file` can either be a genome in the FASTA format or a pre-built
index generated by
```sh
miniprot -t8 -d ref.mpi ref.fna
```
Because miniprot indexing is slow and memory intensive, it is recommended to
pre-build the index. FASTA input files can be optionally compressed with gzip.

Miniprot outputs alignment in the protein PAF format. Different from the more
common nucleotide PAF format, miniprot uses more CIGAR operators to encode
introns and frameshifts. Please refer to the manpage for detailed explanation.
You may also convert the output PAF to GFF3 with
```sh
paftools.js paf2gff -a aln.paf > out.gff
```
This requires the latest paftools.js from the minimap2 github.

### Preliminary evaluation

We aligned 21,919 canonical mouse proteins from Gencode M30 against pre-indexed
GRCh38. It took about 5 minutes over 16 threads. Miniprot mapped 19,309
proteins with 154,660 introns in their best alignment. In comparison to Gencode
v40 human annotations, 94.5% of these introns are annotated in Gencode with
exact coordinates.

We also aligned 52,089 zebrafish proteins, including alternative isoforms,
against GRCh38. Miniprot mapped 34,283 of them with 84.5% of 267,776 predicted
introns confirmed by the Gencode annotation. The accuracy is lower because
zebrafish is more distant from human and because the input includes rarer
isoforms.

## Limitations

* Miniprot may report duplicated alignments for genes in tandem segmental
  duplications.

* The initial conditions of dynamic programming are not technically correct,
  which may result in suboptimal residue alignment in rare cases.

* Support for non-splicing alignment needs to be improved.

* More manual inspection required for improved accuracy.
