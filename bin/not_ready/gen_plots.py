#!/usr/bin/env python3

import scanpy as sc
from argparse import ArgumentParser
from pathlib import Path

# Initialize argument parser, add arguments, and parse
parser = ArgumentParser()
parser.add_argument("--anndata_path", "-a", required=True, type=Path)
args = parser.parse_args()

# Read anndata object into memory
adata = sc.read_h5ad(args.anndata_path)

# Establish plots directory and set the figure directory in sc settings
plots_dir = Path("plots")
plots_dir.mkdir()
sc.settings.figdir = plots_dir

# Define a mapping between different plot types and gene types
plot_colors = {"total_counts": ["hemoglobin", "sex_linked", "stress_response"],
               "pct_counts": ["mitochondrial", "ribosomal", "cell_cycle"]}

# Generate the plots and save
for plot_type, gene_types in plot_colors.items():
    for gene_type in gene_types:
        color = f"{plot_type}_{gene_type}"
        save_path = f"{color}.svg"
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color=color, save=save_path) #type: ignore