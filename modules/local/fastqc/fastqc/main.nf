#!/usr/bin/env nextflow

/*
 * Quality control on FASTQ files
 */

 process FASTQC {

    publishDir params.outdir + '/fastqc', mode: 'copy'

    input:
        path(reads)

    output:
        path '*'    // *.html, *.zip might be better

    script:
    """
    fastqc -t 16 ${reads} 
    """
}
