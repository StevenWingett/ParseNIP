#!/usr/bin/env nextflow

/*
 * QC on the split-pipe mapping results
 */

 process SPLITPIPE_QC {

    publishDir params.outdir + '/scanpy_mapping_QC', mode: 'copy'

    input:
        path(results_folders)

    output:
        path 'mapping_summary_plots'

    script:
    """
    python3 ${projectDir}/bin/parse_mapping_qc.py
    """
}