#!/usr/bin/env nextflow

/*
 * Create a Seruat object using the split-pipe output files
 */

 process SEURAT_BUILDSEURATOBJECT {

    publishDir params.outdir + '/seurat_data_output', mode: 'copy'

    input:
        path(results_folders)
        val(prefix)

    output:
        //path 'seurat_data_output'
        path 'seurat_data_objects_*'

    script:
    """
    R --version
    #ls *
    #mkdir seurat_data_output
    for D in */all-sample/DGE_filtered/; do echo \$D; Rscript ${projectDir}/bin/create_seurat_object.R \$D \$D ${prefix}; done
    """
}