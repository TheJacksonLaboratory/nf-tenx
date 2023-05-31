#!/usr/bin/env python3

from glob import glob
from pathlib import Path

import scanpy as sc
from numpy import loadtxt


def gen_plots() -> None:
    # The only thing staged in by Nextflow should be the anndata file
    # (extension .h5ad). glob for it and read
    adata_path = glob('**/*.h5ad', recursive=True)[0]
    adata = sc.read_h5ad(adata_path)

    # Establish plots directory and set the figure directory in sc
    # settings
    plots_dir = Path('plots')
    plots_dir.mkdir()
    sc.settings.figdir = plots_dir

    # Define a default mapping between plot types and gene types
    plot_colors = {
        'total_counts': {'hemoglobin', 'sex_linked', 'stress_response'},
        'pct_counts': {'mitochondrial', 'ribosomal', 'cell_cycle'},
    }

    # Check that the annotated types are the default, if not just plot 
    # total_counts
    annotated_gene_types = {col for col in adata.var.columns if adata.var[col].dtype == bool}
    if annotated_gene_types != plot_colors['total_counts'] | plot_colors['pct_counts']:
        plot_colors = {'total_counts': annotated_gene_types}

    # Generate the plots and save
    for plot_type, gene_types in plot_colors.items():
        for gene_type in gene_types:
            color = f'{plot_type}_{gene_type}'
            save_path = f'{color}.svg'
            sc.pl.scatter(
                adata,
                x='total_counts',
                y='n_genes_by_counts',
                color=color,
                save=save_path,
            )  # type: ignore


if __name__ == '__main__':
    gen_plots()
