#!/usr/bin/env Rscript

library(SoupX)
library(DropletUtils)

# Load data and estimate soup profile
soup_object <- load10X(getwd())

# Estimate level of background contamination
soup_object <- autoEstCont(soup_object)

# Adjust counts accordingly
output <- adjustCounts(soup_object)

# Glob for the features TSV, which containes gene IDs and symbols.
# There should only be one, so just pick the first element
genes_file <- Sys.glob("*filtered*feature*matrix*/*features*.tsv.gz")

# Read the file into a dataframe (note the sep argument and the header
# argument per the format of the 10x file)
gene_info <- read.csv(genes_file, sep = "\t", header = FALSE)

# A CSV mapping each gene to its genome should have been passed in as a
# command-line argument
command_line_args <- commandArgs(trailingOnly = TRUE)
genes_by_genome <- read.csv(command_line_args[2])

output_dir <- command_line_args[3]

# Use the droplet-utils writing function to write the file,
# manually setting the gene IDs, symbols, and  with the dataframes loaded above
write10xCounts(
    file.path(output_dir, "filtered_matrix.h5"),
    output,
    gene.id = gene_info$V1,
    gene.symbol = gene_info$V2,
    gene.type = gene_info$V3,
    genome = genes_by_genome,
    version = "3"
)

write10xCounts(
    file.path(output_dir, "filtered_matrix"),
    output,
    gene.id = gene_info$V1,
    gene.symbol = gene_info$V2,
    gene.type = gene_info$V3,
    genome = genes_by_genome,
    version = "3"
)
