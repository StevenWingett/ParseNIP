#!/usr/bin/env nextflow

/*
 * Generate Parse split-pipe index file
 */
process SPLITPIPE_INDEX {

    publishDir params.outdir, mode: 'copy'

    input:
        path input_fasta
        path input_gtf
        val genome_name

    output:
        path "GENOME_INDEX"

    script:
    """
    split-pipe --version
    split-pipe \
    --mode mkref \
    --genome_name ${genome_name} \
    --fasta ${input_fasta} \
    --genes ${input_gtf} \
    --output_dir GENOME_INDEX
    """
}
