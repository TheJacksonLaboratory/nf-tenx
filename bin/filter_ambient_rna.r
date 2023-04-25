#!/usr/bin/env Rscript

library(SoupX)
library(DropletUtils)

# Load data and estimate soup profile
soup_object <- load10X(getwd())

# Estimate level of background contamination
soup_object <- autoEstCont(soup_object)

# Adjust counts accordingly
output <- adjustCounts(soupx_object)

# Glob for the features TSV, which containes gene IDs and symbols.
# There should only be one, so just pick the first element
genes_file <- Sys.glob("**/*filtered*feature*matrix*/*features*.tsv.gz")[0]
genes_by_genome <- some_file

# Read the file into a dataframe (note the sep argument and the header
# argument per the format of the 10x file)
gene_info <- read.csv(genes_file, sep = "\t", header = FALSE)
genome_by_gene <- read.csv()

# Use the droplet-utils writing function to write the file,
# manually setting the gene IDs/symbols with the dataframe loaded above
write10xCounts("soupx_filtered_matrix.h5",
    output,
    gene.id = gene_ids$V1,
    gene.symbol = gene_ids$V2,
    gene.type = gene_info$V3
)
file.path()