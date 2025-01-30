#!/usr/bin/env nextflow

/*
 * Generate Parse split-pipe index file
 */
process SPLITPIPE_INDEX {

    publishDir params.outdir, mode: 'copy'

    input:
        path input_fasta
        path input_gtf
        val genome_name

    output:
        path "GENOME_INDEX"

    script:
    """
    # Generate parameter file
    # Increase the allowed RAM and the number of junctions to be inserted on the
    # fly to enable the processing of very large genomes (e.g. multi-species)
    echo "mkref_star_args --limitGenomeGenerateRAM,100000000000,--limitSjdbInsertNsj,2000000" > parfile.txt
    
    split-pipe --version
    
    split-pipe \
    --mode mkref \
    --genome_name ${genome_name} \
    --fasta ${input_fasta} \
    --genes ${input_gtf} \
    --output_dir GENOME_INDEX \
    --parfile parfile.txt
    """
}
