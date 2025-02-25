#!/usr/bin/env nextflow

/*
 * Quality control on FASTQ files
 */

 process FASTQC {

    publishDir params.outdir + '/fastqc', mode: 'copy'

    input:
        path(reads)
        val control   // Used to prevent this process executing until checks are finished

    output:
        path '*'    // *.html, *.zip might be better

    script:
    """
    fastqc -t 16 ${reads} 
    """
}
