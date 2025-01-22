#!/usr/bin/env nextflow

/*
 * MultiQC summary of results
 */

 process MULTIQC {

    publishDir params.outdir + '/multiqc', mode: 'copy'

    input:
        path(results_folders)

    output:
        path '*.html'

    script:
    """
    multiqc *
    """
}
