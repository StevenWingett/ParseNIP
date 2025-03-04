#!/usr/bin/env nextflow

process CHECK_SETUP_PREBUILT_INDEX {

    input:
        path(fastq_folder)
        path(samp_list)
        val(chemistry)
        path(genome_dir)

    output:
        val('CONTROL_1')    // Prevents later processes executing until after the check is completed

    script:
    """
    python3 ${projectDir}/bin/check_setup.py ${fastq_folder} ${samp_list} ${chemistry} --genome_dir ${genome_dir}
    """
}