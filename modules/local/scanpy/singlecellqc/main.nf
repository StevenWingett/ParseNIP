#!/usr/bin/env nextflow

/*
 * Scanpy Single Cell QC on filtered matrix file
 */

 process SCANPY_SINGLECELLQC {

    publishDir params.outdir + '/scanpy_SingleCellQC', mode: 'copy'

    input:
        path(results_folders)
        val(prefix)

    output:
        path 'scanpy_singlecellqc_*'


    script:
    """
    for D in */*/DGE_filtered/; do echo \$D; python3 ${projectDir}/bin/scanpy_qc.py --directory \$D --outdir \$D --pipeline ${prefix}; done
    """
}