#!/usr/bin/env Rscript

library(Seurat)
library(SeuratObject)

# Get the conversion directory from the command line
command_line_args <- commandArgs(trailingOnly = TRUE)
total_counts_dir <- command_line_args[2]
output_path <- paste(total_counts_dir, "_annotated.rds", sep = "")

# Read in raw data and gene/cell metadata
raw_data <- Read10X(total_counts_dir)
obs <- read.csv("obs.csv")
var <- read.csv("var.csv")


# Create the Seurat object from raw data and assign the gene/cell annotations
# appropriately
sobj <- CreateSeuratObject(counts = raw_data, meta.data = obs)
sobj@assays$RNA.var.features <- var

dirs <- list.dirs(recursive = FALSE)
splice_dirs <- dirs[file.path(dirs) != total_counts_dir]
for (dir in splice_dirs) {
    # Read the raw data in and create an assay object from it
    data <- Read10X(data.dir = dir)
    aobj <- CreateAssayObject(counts = data)

    # Set a new assay in the existing Seurat object
    sobj[[dir]] <- aobj
}

saveRDS(sobj, file = output_path)
