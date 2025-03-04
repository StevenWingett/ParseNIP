#!/usr/bin/env nextflow

process CHECK_SETUP_PREBUILT_INDEX {

    input:
        //val(fastq_folder_name)
        path(fastq_folder)
        //val(samp_list_name)
        path(samp_list)
        val(chemistry)
        //val(genome_dir_name)
        path(genome_dir)

    output:
        val('CONTROL_1')    // Prevents later processes executing until after the check is completed

    script:
    """
    python3 ${projectDir}/bin/check_setup.py ${fastq_folder} ${samp_list} ${chemistry} --genome_dir ${genome_dir}
    """
}