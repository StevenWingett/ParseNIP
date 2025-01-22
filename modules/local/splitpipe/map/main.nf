#!/usr/bin/env nextflow

/*
 * Map FASTQ file data with split-pipe
 */

process SPLITPIPE_MAP {

    publishDir params.outdir + "/splitpipe", mode: 'copy'

    input:
        tuple val(sublibrary_id), path(reads)
        path genome_dir
        path samp_list
        val chemistry

    output:
        path sublibrary_id
    
    script:    // add --dryrun for testing purposes
    """
    echo "post_min_map_frac 0.001" > parfile.txt

    split-pipe \
    --mode all \
    --chemistry ${chemistry} \
    --genome_dir ${genome_dir} \
    --fq1 ${reads[0]} \
    --fq2 ${reads[1]} \
    --output_dir ${sublibrary_id} \
    --parfile parfile.txt \
    --samp_list ${samp_list} \
    """
}