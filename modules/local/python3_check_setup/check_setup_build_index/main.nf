#!/usr/bin/env nextflow

process CHECK_SETUP_BUILD_INDEX {

   input:
        path(fastq_folder)
        path(samp_list)
        val(chemistry)
        path(fasta_files)
        path(gtf_files)
        val(genome_names)

    output:
        val('CONTROL_1')    // Prevents later processes executing until after the check is completed

    script:
    """   
    python3 ${projectDir}/bin/check_setup.py ${fastq_folder} ${samp_list} ${chemistry} --fasta_files ${fasta_files} --gtf_files ${gtf_files} --genome_names ${genome_names} --basename
    """
}