#!/usr/bin/env nextflow

/*
 * Combine the separate split-pipe sublibraries results
 */

 process SPLITPIPE_COMBINE {
    
    publishDir params.outdir, mode: 'copy'

    input:
        path sample_ids
        path genome_dir     // The genome path is needed to make these files visible
                            // (the path is actually stored by the split-pipe mapping output files)
    output:
        path "splitpipe_combined_sublibraries"

    script:   // Add --dryrun for testing purposes
    """
    split-pipe --mode comb --sublibraries ${sample_ids} --output_dir splitpipe_combined_sublibraries
    """
}