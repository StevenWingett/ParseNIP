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
    # MultiQC option reports 2 directories to avoid name clashes with barcode_headAligned_anno.bam
    multiqc  -d -dd 2 *
    """
}
