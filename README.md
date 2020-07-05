README
======

ProDA (Protein Domain Aligner) is public
domain software for generating multiple alignments of protein sequences
with repeats and rearrangements, e.g. proteins with multiple domains.

Given a set of protein sequences as
input, ProDA first finds local pairwise alignments between all pairs of
sequences, then forms blocks of alignable sequence fragments, and
finally generates multiple alignments of the blocks. ProDA relies on
many techniques used in Probcons
([http://probcons.stanford.edu](http://probcons.stanford.edu/)), a
recent multiple aligner that shows high
accuracy in a number of popular benchmarks.

Algorithm developed by Tu Minh Phuong,
Chuong B. Do, Robert C. Edgar, and Serafim Batzoglou.


# Algorithm outline

1. Compute local pairwise alignments for each pair of
   sequences using either Viterbi or posterior decoding.
2. Infer repeats from pairwise alignments.
3. Generate a block of *alignable sequence fragments.*
4. Construct guide tree using *expected accuracies* and adjustment of the block.
5. *Progressively align* the block using the guide tree.
6. Extract final alignments from block alignment.
7. Remove used pairwise alignments.


# Usage

To install and use ProDA, `make` and use as follows

```bash
./proda input.fasta > output.aln
```

Provide the input via in FASTA format. Output will be written
in Clustal format.

For a set of input sequences, Proda 
usually outputs several blocks in turn, each consists of alignable
sequence fragments. Each block is followed by its multiple
alignment. 

A block is specified by listing its sequence fragments.
Each fragment is output as `sequence_name(start-end)`.


# Command line options

* `-L min_length` - Set minimal alignment length equal to min_length.

ProDA finds alignments of length greater than or equal to a
threshold *LMIN*. By default, *LMIN* = 30. This
option sets the threshold to min_length.

* `-posterior` - Use posterior decoding when computing local pairwise alignments

ProDA computes local pairwise alignments between two sequences
using a pair-HMM and either Viterbi decoding or posterior decoding. The
default option is using Viterbi decoding which is faster than posterior
decoding but may be less accurate. Turning on this option instructs the
aligner to use posterior decoding instead.

* `-silent` - Do not report progress while aligning.

Turning on this option instructs the aligner not to report the
progress while aligning. By default, ProDA reports the progress on all
pairwise alignments, block generation, and on block alignment. Since
some stages of the algorithm, especially pairwise alignment, may take
long time, reporting progress makes the program look alive while
running.

* `-tran` - Use transitivity when forming blocks of alignable sequence fragments.

Two sequence fragments are *directly* *alignable* if they are
parts of a local pairwise alignment. By default, two fragments are
considered *alignable* if and only if they are directly alignable.
Turning on this option instructs the aligner to consider two fragments
alignable when they are directly alignable or when both of them are
directly alignable to a third fragment.

* `-fasta` - Use FASTA output format in addition to the ClustalW format.

When this option is turned on, the
aligner generates output in the FASTA format and stores in a file with
the same name as the first input file and extension “.fasta”, in
addition to the normal output to stdout.
Aligned residues are in upper case, unaligned residues are in lower case.
This option should be used only if input sequences do not contain repeats.

