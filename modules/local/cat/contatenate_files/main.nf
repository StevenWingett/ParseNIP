#!/usr/bin/env nextflow

process CONCATENATE_FILES {

    publishDir params.outdir + "/concatenated_files", mode: 'copy'

    input:
        tuple val(id), path(f_fastqs), path(r_fastqs)

    output:
        tuple val(id), path("*.fq.gz")

    script:
    """
    cat $f_fastqs > ${id}.r_1.fq.gz
    cat $r_fastqs > ${id}.r_2.fq.gz
    """
}