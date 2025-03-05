# Parse Nextflow-Integrated Pipeline (ParseNIP) <br> *(Processing Parse Evercode data using Split-Pipe)* <br> ![GitHub release (latest by date)](https://img.shields.io/github/v/release/StevenWingett/ParseNIP)

## Overview
This pipeline is designed for the processing of Evercode Datasets using split-pipe as the aligner.  Split-pipe was developed by [Parse Biosciences](https://www.parsebiosciences.com).

Only 1 Parse "experiment" should be processed at a time.  Each experiment can include replicates, sublibraries etc.  There should only be 1 genome reference (although this may comprise multiple species.)  The pipeline expects **paired-end** sequence data.

Pipeline steps:

1. Build a genome index (optional: `--fasta` `--gtf` `--genome_name`).  To do this, specify a FASTA file, a GTF file and a genome name.  Alternatively, the FASTA and GTF files may be specified using the `--genome` option (see: https://nf-co.re/docs/usage/reference_genomes).  The `--genome` option uses information in a nextflow configuration file to determine the correct FASTQ and GTF files to use.

    If you do NOT wish to create a genome de-novo, but would rather use a pre-existing genome, then specify the `--genome_dir` option.  REMEMBER: Parse recommends that a genome index should be built by the same version of split-pipe that is used to perform the mapping.  Thus, building the genome de-novo each time the pipeline is run ensures that this is the case.

    If you need to build a multi-species genome, please specify the FASTA files, GTF files and genome names as a comma-separated list.  For example, to build an index comprising 2 genomes:

    `--fasta genome1.fasta,genome2.fasta --gtf genome1.gtf,genome2.gtf --genome_name genome1_name,genome2_name`

    It is also possible to instruct the pipeline to only build a genome index i.e. to NOT process FASTQ files.  To do that, specify the parameters to build the index (`--fasta` /`--gtf` / `--genome_name`), but do not specify other parameters (e.g. `--fastq`). 

2.  FASTQC - General quality check on the sequence data.  (This step can be omitted with the option `--skip_fastqc`.)

3.  Trim Galore (optional: `--trim`) - hard-trim input sequence reads to 64bp.  Very long reads (e.g. >100 bp are more likely to contain non-transcriptomic sequence at the 3' end).  [In future updates the trimming feature may be extended to quality trim the sequences.]

4.  Concatenate (optional: `--concatenate`) - FASTQ files from the sample sublibrary (sometimes these datasets may be split across multiple lanes).

5.  Map the reads using the mapper split-pipe, produced by Parse Biosciences.  To build the genome de-novo: `--fasta` `--gtf` `--genome_nam`e; to reference a pre-build folder: `genome_dir`; to reference a pre-build genome using a configuration file: `--genome` `--genome_name`. 

6.  Combine the results from different sublibraries (if present).

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

Here we have 2 sub-libraries, "test" and "test2".  Sublibrary "test" was run on lane1, whereas sublibrary2 was run on lane 1 and lane 2.  The option --concatenate will join the FASTQ files for test2 prior to mapping.  So, sublibrary files are merged before mapping, but different sublibraries are kept separate.  After mapping however, separate sublibraries are combined using a split-pipe algorithm.

## A note on the sample list file
This should be a **space-delimited file**, listing each sample and well ID e.g.

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

This is the format that split-pipe takes as input if being run outside of nextflow.

## Example command
The following command illustrates how to use the nf-split-pipe nextflow pipeline.  The processing has been configured using the file nextflow.config, using a profile named lmb_cluster.  The space-delimited file parse_samplesheet.txt details the plate well / sample relationships.  Prior to mapping, the FASTQ reads will be trimmed and FASTQ files from the same sublibrary will be concatenated.  The job will be run in the background (`-bg`) - we strongly advise running nf-split-pipe jobs in the background, as processing may take hours/days to complete.  

`nextflow run -config nextflow.config -profile lmb_cluster ./nf-split-pipe/main.nf --samp_list parse_samplesheet.txt --fastq FASTQ --trim --concatenate --gtf Homo_sapiens.GRCh38.102.gtf --fasta homo_sapiens__GRCh38__release102.dna.fa --genome_name Homo_sapiens.GRCh38.102 -bg`

## GUIde-Piper
Remember that pipeline commands can be built using [GUIde-Piper](http://guidepiper/parse), within the LMB.

2025 Steven Wingett, The MRC-Laboratory of Molecular Biology, Cambridge, UK.
