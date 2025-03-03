#!/usr/bin/env nextflow

process CHECK_SETUP_BUILD_INDEX {

   input:
        //val(fastq_folder_name)
        path(fastq_folder)
        //val(samp_list_name)
        path(samp_list)
        val(chemistry)
        //val(fasta_files_names)
        path(fasta_files)
        //val(gtf_files_names)
        path(gtf_files)
        val(genome_names)

    output:
        val('CONTROL_1')    // Prevents later processes executing until after the check is completed

    //# python3 ${projectDir}/bin/check_setup.py ${fastq_folder_name} ${samp_list_name} ${chemistry} --fasta_files ${fasta_files_names} --gtf_files ${gtf_files_names} --genome_names ${genome_names} --basename
    script:
    """
    echo test > output.txt

    
    python3 ${projectDir}/bin/check_setup.py ${fastq_folder} ${samp_list} ${chemistry} --fasta_files ${fasta_files} --gtf_files ${gtf_files} --genome_names ${genome_names} --basename
    """
}