#!/usr/bin/env nextflow

process CHECK_SETUP_PREBUILT_INDEX {

    publishDir params.outdir + "/test", mode: 'copy'

    //input:
    //    val(id)

    output:
        val('CONTROL_1')    // Prevents later processes executing until after the check is completed


    script:
    """
    sleep 10
    """
}