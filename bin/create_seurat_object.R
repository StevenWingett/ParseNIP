# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # Rscript based an the R script produce by Parse Biosciences: 
# # https://support.parsebiosciences.com/hc/en-us/articles/360053078092-Seurat-Tutorial-65k-PBMCs
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



#rm(list = ls())


# Setup

# directory: Path to the split-pipe output directory that needs processing
# outdir: Write result to here
# pipeline: Script being used as part of pipeline - if so use custom output directory location
args = commandArgs(trailingOnly=TRUE)

# Check 3 arguments passed
if (length(args)!=3) {
  stop("Error: provide 3 arguments: directory, outdir, pipeline\n", call.=FALSE)
}


#data_path <- "/mnt/seurat_data_output/"
#mat_path <- "/mnt/dan_test_data"
data_path <- args[2]   # Output directory
data_path <- strsplit(data_path, "/")
data_path <- data_path[[1]]
data_path <- data_path[1:length(data_path) - 1]
data_path <- paste(data_path, collapse="/")
data_path <- paste0('seurat_data_objects_', args[3], '_', data_path, '/')

mat_path <- args[1]    # Input directory


# Convenience functions
SaveObject <- function(object, name){
  saveRDS(x, paste0(data_path, name, ".RDS"))
}

x <- c(1, 2, 3)


dir.create(data_path, recursive=TRUE)

SaveObject( x, "seurat_obj_before_QC")
#ls 

print(paste("input:", mat_path))
print("")
print(paste("output:", data_path))



