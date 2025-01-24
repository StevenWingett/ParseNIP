# Parse Evercode Split-Pipe Nextflow Pipeline

## Overview
This pipeline is designed for the processing of Evercode Datasets using split-pipe as the aligner.  Split pipe is developed by [Parse Biosciences](https://www.parsebiosciences.com).

Only 1 Parse "experiment" should be processed at a time.  Each experiemnt can include replicates, sublibraries etc.  The pipeline expects **paired-end** sequence data.

Pipeline steps:

1.  FASTQC - General quality check on the sequence data

2.  Trim Galore (optional: --trim) - hardtrim input sequence reads to 64bp.  Very long reads (e.g. >100 bp are more likely to contain non-transriptomic sequence at the 3' end).

3.  Concatenate (optional: --concatenate) - FASTQ files from the sample sublibrary (sometimes these datasets may be split across multiple lanes).

4.  Build a genome index (optional).  To do this, specify a --fasta file, a --gtf file and a --genome_name.  Alternatively, the FASTA file and GTF file may be specified using the --genome option, which takes a genome config file (see: https://nf-co.re/docs/usage/reference_genomes).  If you do NOT wish to create a genome denovo, but would rather use a pre-existing genome, then specify the --genome_dir option.  REMEMBER: Parse recommends that a genome index should be built by the same version of split-pipe that is used to perform the mapping.  Thus, building the genome denovo each time the pipeline is run ensures that this is the case.

If you need to build a multi-species genome, please specify the FASTA files, GTF files and genome names as a comma-separted list.  For example, to build an index comprising 2 genomes:

`--fasta genome1.fasta,genome2.fasta --gtf genome1.gtf,genome2.gtf --genome_name genome1_name,genome2_name`

5.  Map the reads using the mapper split-pipe, produced by Parse Biosciences.

6.  Combine the results from different sub-libraries (if present).

7.  Summarise results in a MultiQC report

## A note on file naming

Filenames need to be of the format:
[identifier].s_[lane number].r_[read1 or read2].fq.gz

For example: \
test1.s_1.r_1.fq.gz \
test1.s_1.r_2.fq.gz \
test2.s_1.r_1.fq.gz \
test2.s_1.r_2.fq.gz \
test2.s_2.r_1.fq.gz \
test2.s_2.r_2.fq.gz 

Here we have 2 sub-libraries, "test" and "test2".  Sublibrary "test" was run on lane1, whereas sublibrary2 was run on lane 1 and lane 2.  The option --concatenate will join (i.e. concatenate!) the FASTQ files for test2 prior to mapping.  So, sublibrary files are merged before mapping, but different sublibraries are kept separate.  After mapping however, separate sublibraries are combined using a split-pipe algorithm.

## A note on the sample list file
This should be a space delimited file, listing each sample and well ID e.g.

Sample1 A1 \
Sample2 A2 \
Sample3 A3 \
Sample4 A4 \
Sample5 A5 \
Sample6 A6 \
Sample7 A7 \
Sample8 A8 \
Sample9 A9 \
Sample10 A10 \
Sample11 A11 \
Sample12 A12 \
Sample13 B1 \
Sample14 B2 \
. \
. \
. 