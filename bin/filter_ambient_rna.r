#!/usr/bin/env Rscript

library(SoupX)
library(DropletUtils)
library(rhdf5)

# Get command-line arguments and store in variables
command_line_args <- commandArgs(trailingOnly = TRUE)
cellranger_version <- command_line_args[2]
cellranger_dir <- command_line_args[3]
soupx_dir <- command_line_args[4]

# Load data, estimate soup profile and background contamination, adjust counts
soup <- load10X(cellranger_dir)
soup <- autoEstCont(soup)
filtered_counts <- adjustCounts(soup)

# Read the features.tsv file into a dataframe
genes_file <- Sys.glob("**/*filtered*feature*matrix*/*features*.tsv.gz")[1]
gene_info <- read.csv(genes_file, sep = "\t", header = FALSE)

# Get first digit of cell ranger version
cellranger_version <- strsplit(
    cellranger_version,
    split = ".",
    fixed = TRUE
)[[1]][1]

# write10XCounts takes in either "2" or "3" for the version parameter, so set
# cellranger_version manually because if it is not "2" or "3", write10XCounts
# just defaults to "2"
if (cellranger_version >= "3") {
    cellranger_version <- "3"
} else {
    cellranger_version <- "2"
}

# Get the genomes out of the molecule_info.h5 file that cellranger generates
mol_file <- Sys.glob("**/*molecule*info*.h5")[1]
mol_info <- H5Fopen(mol_file)
genomes <- mol_info$features$genome

# Cellranger's filtered matrix comes as a .h5 file and as a matrix market
# directory. Loop over the corresponding extensions and save the soupx matrix
# with the same name as cellranger's
extensions <- c("*.h5", "")
for (ext in extensions) {
    # Generate the glob pattern and glob
    pattern <- paste("**/*filtered*matrix", ext, sep = "")
    cr_filt_mtx <- Sys.glob(pattern)[1]

    # Get the name of the cellranger file and save the soupx matrix to the
    # soupx_dir with the same name
    file_name <- basename(cr_filt_mtx)
    soupx_path <- file.path(soupx_dir, file_name)
    write10xCounts(
        soupx_path,
        filtered_counts,
        gene.id = gene_info$V1,
        gene.symbol = gene_info$V2,
        gene.type = gene_info$V3,
        genome = genomes,
        version = cellranger_version
    )
}
H5Fclose(mol_info)