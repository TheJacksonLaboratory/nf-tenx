#!/usr/bin/env python3

#random change
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import scanpy as sc
import scrublet as scr

# Command line arguments should be pickled dataframe paths.
# The dataframes contain genes, annotated by type
parser = ArgumentParser()
parser.add_argument('--ref_df_paths', '-r', nargs='*', type=Path, required=True)
args = parser.parse_args()

# The current directory should have cellranger outputs.
# Get the filtered feature matrix, ensuring there is just one.
filt_mat_paths = list(Path('.').glob('*filtered*matrix*.h5'))
if len(filt_mat_paths) != 1:
    raise ValueError(
        f'The following filtered matrix files were found:\n{", ".join(str(path) for path in filt_mat_paths)}\nThis pipeline requires exactly one.'
    )
adata = sc.read_10x_h5(filt_mat_paths[0])

# Iterate over each reference dataframe and each gene type and annotate
for df_path in args.ref_df_paths:
    ref_df = pd.read_pickle(df_path)

    # The boolean columns of the reference dataframes are the gene
    # annotations, so get them with generator comprehension
    gene_types = (col for col in ref_df.columns if ref_df[col].dtype == bool)
    
    # Iterate over each gene type, creating a column in adata for each
    # that stores a boolean for whether the gene is of that type or not
    for gene_type in gene_types:
        adata.var[gene_type] = (
            adata.var['gene_ids']  # Match by gene ID
            .str.upper()  # Convert gene IDs to uppercase
            .isin(ref_df.loc[ref_df[gene_type], 'ensembl_gene_id'])
        )  # If the gene ID matches those in the gene type column of
           # the reference dataframe, mark the gene as true for that
           # gene type

# Make genes/cells names unique
adata.var_names_make_unique()
adata.obs_names_make_unique()

# Calculate QC metrics. Generate a list of annotated gene types first
annotated_gene_types = [
    col for col in adata.var.columns if adata.var[col].dtype == bool
]
sc.pp.calculate_qc_metrics(adata=adata, qc_vars=annotated_gene_types, inplace=True)

# Calculate predicted doublets and doublet scores
scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata.obs['doublet_score'] = doublet_scores
adata.obs['doublet_predicted'] = predicted_doublets

# Save to anndata object to disk
output_path = Path('final_anndata.h5ad')
adata.write_h5ad(output_path)
