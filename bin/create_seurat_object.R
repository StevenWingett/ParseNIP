# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # Rscript based an the R script produce by Parse Biosciences: 
# # https://support.parsebiosciences.com/hc/en-us/articles/360053078092-Seurat-Tutorial-65k-PBMCs
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Setup
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)


# directory: Path to the split-pipe output directory that needs processing
# outdir: Write result to here
# pipeline: Script being used as part of pipeline - if so use custom output directory location
args = commandArgs(trailingOnly=TRUE)

# Check 3 arguments passed
if (length(args)!=3) {
  stop("Error: provide 3 arguments: directory, outdir, pipeline\n", call.=FALSE)
}

data_path <- args[2]   # Output directory
data_path <- strsplit(data_path, "/")
data_path <- data_path[[1]]
data_path <- data_path[1:length(data_path) - 1]
data_path <- paste(data_path, collapse="/")
data_path <- paste0('seurat_data_objects_', args[3], '_', data_path, '/')

mat_path <- args[1]    # Input directory


# Convenience functions
SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

mat <- ReadParseBio(mat_path)


# Check to see if empty gene names are present, add name if so.
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"


# Read in cell meta data
cell_meta <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)


# Create object
resultsObject <- CreateSeuratObject(mat, min.features = 0, min.cells = 0,   # In the original script from Parse, min.features is min.genes (and threshold was set > 0)
names.field = 0, meta.data = cell_meta)


# Setting our initial cell class to a single type, this will changer after clustering. 
resultsObject@meta.data$orig.ident <- factor(rep("resultsObject", nrow(resultsObject@meta.data)))
Idents(resultsObject) <- resultsObject@meta.data$orig.ident


dir.create(data_path, recursive=TRUE)  # Create output directory   
SaveObject(resultsObject, "seurat_obj_before_QC")

print('Done')
