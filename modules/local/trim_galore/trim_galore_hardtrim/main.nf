#!/usr/bin/env nextflow

/*
 * Trim FASTQ files
 */


process TRIMGALORE_TRIM {

    publishDir params.outdir + '/trim_galore', mode: 'copy'

    input:
        tuple val(read_id), path(reads)
        val control   // Used to prevent this process executing until checks are finished
        

    output:
        tuple val(read_id), path("*.r_1.64bp_5prime.fq.gz"), path("*.r_2.64bp_5prime.fq.gz"), emit: trimmed_seqs
        // could add an output channel for fastqc here, as it doesn't need to be in the same channel as the reads. 
    
    script:
    """
    trim_galore \
    --cores 4 \
    --paired \
    --hardtrim5 64 \
    ${reads[0]} \
    ${reads[1]} \
    """
}
