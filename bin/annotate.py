#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
import scanpy as sc
import scrublet as scr

# Define the file extensions for each matrix type
main_mat_h5 = '.h5'
main_mat_h5ad = '.h5ad'
velo = '.loom'
ref_data = '.pickle'

# Define a dictionary of readers for each extension and a dictionary to store the items read
readers = {main_mat_h5: sc.read_10x_h5,
           main_mat_h5ad: sc.read_h5ad,
           velo: sc.read_loom,
           ref_data: pd.read_pickle}

# Read files into a dict of lists using dict-comprehension
mat_by_ext = {
    ext: [reader(path) for path in Path.cwd().glob(f'**/*{ext}')]
    for ext, reader in readers.items()
}

main_mats = mat_by_ext[main_mat_h5] + mat_by_ext[main_mat_h5ad]
if len(main_mats) != 1:
    raise ValueError(
        f'Something went wrong. The following filtered feature matrix files were found, but at this point of the '
        f'pipeline, there should just be one.\n'
        + '\n'.join(str(adata.filename) for adata in main_mats))

adata = main_mats[0]

# Add velocyto layers to adata if it's there
if mat_by_ext[velo]:
    velo_adata = mat_by_ext[velo][0]
    for l in (layer for layer in velo_adata.layers if layer != 'matrix'):
        adata.layers[l] = velo_adata.layers[l].copy()

# Iterate over each reference dataframe and each gene type and annotate
for ref_df in mat_by_ext[ref_data]:

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
sc.pp.calculate_qc_metrics(
    adata=adata, qc_vars=annotated_gene_types, inplace=True)

# Calculate predicted doublets and doublet scores
scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata.obs['doublet_score'] = doublet_scores
adata.obs['doublet_predicted'] = predicted_doublets

# Save to anndata object to disk
output_path = Path('final_anndata.h5ad')
adata.write_h5ad(output_path)
