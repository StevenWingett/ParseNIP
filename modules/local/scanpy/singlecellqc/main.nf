#!/usr/bin/env nextflow

/*
 * Scanpy Single Cell QC on filtered matrix file
 */

 process SCANPY_SINGLECELLQC {

    publishDir params.outdir + '/scanpy_SingleCellQC', mode: 'copy'

    input:
        path(results_folders)

    output:
        path '*.txt'
        path 'scanpy_singlecellqc_*'


    script:
    """
    ls -d */*/DGE_filtered/ > DGE_folders_to_process.txt
    for D in */*/DGE_filtered/; do echo \$D; python3 ${projectDir}/bin/scanpy_qc.py --directory \$D --outdir \$D --pipeline; done
    """
}