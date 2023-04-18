#!/usr/bin/env Rscript

library(Seurat)
library(SeuratObject)
library(glue)

# Get the conversion directory from the command line
command_line_args <- commandArgs(trailingOnly=TRUE)

rna_velo <- command_line_args[2]
conversion_dir <- Sys.glob('*')[1]

# Read in raw data and gene/cell metadata
raw_data <- Read10X(data.dir=glue('{conversion_dir}/total_counts'))
obs <- read.csv(glue('{conversion_dir}/obs.csv'))
var <- read.csv(glue('{conversion_dir}/var.csv'))


# Create the Seurat object from raw data and assign the gene/cell annotations appropriately
sobj <- CreateSeuratObject(counts=raw_data, meta.data=obs)
sobj@assays$RNA.var.features <- var

# Save object if RNA velocity not provided
if (is.na(rna_velo))
{
    saveRDS(sobj, file="final_seurat_object")
} else # If it is provided, iterate over splice types and add their matrices as assays to the object
{
    # Define splice types in list and iterate over them
    splice_types <- list('spliced', 'unspliced', 'ambiguously_spliced')
    
    for (splice_type in splice_types)
    {
        # Read the raw data in and create an assay object from it
        data <- Read10X(data.dir = glue('{conversion_dir}/{splice_type}'))
        aobj <- CreateAssayObject(counts=data)
        
        # Set a new assay in the existing Seurat object
        sobj[[splice_type]] <- aobj
    }
    
    saveRDS(sobj, file="final_seurat_object")
}