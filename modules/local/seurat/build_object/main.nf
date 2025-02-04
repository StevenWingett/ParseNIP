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
        path 'seurat_data_objects_*'

    script:
    """
    R --version

    # Only process matrix if fewer than 2 billion observations
    maximum_allowed=2000000000

    for D in */all-sample/DGE_filtered/
    do
      echo \$D;

      observations=\$(head -3 \$D/count_matrix.mtx | tail -1 | cut -d' ' -f3)

      if [ "\$observations" -le "\$maximum_allowed" ]; then
        Rscript ${projectDir}/bin/create_seurat_object.R \$D \$D ${prefix};
      else
        mkdir -p seurat_data_objects_${prefix}_\$D/
        echo "Matrix too large to generate a Seurat object!" > seurat_data_objects_${prefix}_\$D/README.txt;
      fi
    done
    """
}