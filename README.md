# Parse Evercode Split-Pipe Nextflow Pipeline

## Overview
This pipeline is designed for the processing of Evercode Datasets using split-pipe as the alinger.  Split pipe is developed by [Parse Biosciences](https://www.parsebiosciences.com).

Only 1 Parse "experiment" should be processed at a time.  Each experiemnt can include replicates, sublibraries etc.  The pipeline expects **paired-end** sequence data.

Pipeline steps:

1.  FASTQC - General quality check on the sequence data
2.  Trim Galore (optional: --trim) - hardtrim input sequence reads to 64bp.  Very long reads (e.g. >100 bp are more likely to contain non-transriptomic sequence at the 3' end).
3.  Concatenate (optional: --concatenate) - FASTQ files from the sample sublibrary (sometimes these datasets may be split across multiple lanes)
4.  Map the reads using the mapper split-pipe, produced by Parse Biosciences
5.  Combine the results from different sub-libraries (if present)
6.  Summarise results in a MultiQC report

## A note on file naming

Filenames need to be of the format:
[identifier].s_[lane number].r_[read1 or read2].fq.gz

For example:
test.s_1.r_1.fq.gz
test.s_1.r_2.fq.gz
test2.s_1.r_1.fq.gz
test2.s_1.r_2.fq.gz
test2.s_2.r_1.fq.gz
test2.s_2.r_2.fq.gz

Here we have 2 sub-libraries, "test" and "test2".  Sublibrary "test" was run on lane1, whereas sublibrary2 was run on lane 1 and lane 2.  The option --concatenate will join (i.e. concatenate!) the FASTQ files for test2 prior to mapping.  So, sublibrary files are merged before mapping, but different sublibraries are kept separate.  After mapping however, separate sublibraries are combined using a split-pipe algorithm.


