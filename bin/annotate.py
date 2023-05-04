#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
import scanpy as sc
import scrublet as scr
from arg_utils import parse_cl_paths


def annotate(output_path: Path) -> str:
    """
    This function takes a filtered feature matrix file and annotates it,
    saving the result in the input parameter output_path.

    It annotates various gene types as well as doublet scoring by scrublet,
    also calculating QC metrics while it's at it.
    """
    # Define glob patterns for each matrix type
    filt = '*filtered*matrix*.h5'
    velo = '*.loom'
    ref = '*.parquet'

    # Define a dictionary of readers for each extension and a dictionary
    # to store the items read
    readers = {filt: sc.read_10x_h5, velo: sc.read_loom, ref: pd.read_parquet}

    # Read files into a dict of lists using dict-comprehension
    mtx_by_ext = {
        pattern: tuple(reader(path) for path in Path().rglob(f'{pattern}'))
        for pattern, reader in readers.items()
    }

    too_many = {
        pattern: mtx_list
        for pattern, mtx_list in mtx_by_ext.items()
        if pattern != ref and len(mtx_list) > 1
    }

    if too_many:
        raise ValueError(
            f'Something went wrong. The following is a list of glob patterns and a '
            f'list of files that matched each glob pattern. At this point of the '
            f'pipeline, only one file should match each pattern:\n'
            + '\n'.join(too_many)
        )

    adata = mtx_by_ext[filt][0]

    # Add velocyto layers to adata if it's there
    if mtx_by_ext[velo]:
        velo_adata = mtx_by_ext[velo][0]
        for l in (layer for layer in velo_adata.layers if layer != 'matrix'):
            adata.layers[l] = velo_adata.layers[l].copy()

    # Initialize a set to store all of the gene types in case the
    # reference dataframes have different columns
    all_gene_types = set()

    # Iterate over each reference dataframe and each gene type and
    # annotate
    for ref_df in mtx_by_ext[ref]:
        # The boolean columns of the reference dataframes are the gene
        # annotations, so get them with set comprehension
        current_gene_types = {
            col for col in ref_df.columns if ref_df[col].dtype == bool
        }
        all_gene_types |= current_gene_types  # type: ignore

        # Iterate over each gene type, creating a column in adata for
        # each that stores a boolean for whether the gene is of that
        # type or not
        for gene_type in current_gene_types:
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

    # Calculate QC metrics on annotated gene types
    sc.pp.calculate_qc_metrics(adata=adata, qc_vars=all_gene_types, inplace=True)

    # Calculate predicted doublets and doublet scores
    scrub = scr.Scrublet(adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    adata.obs['doublet_score'] = doublet_scores
    adata.obs['doublet_predicted'] = predicted_doublets

    # Save anndata object to h5ad
    adata.write(output_path)

    return (
        f'The AnnData object associated with {output_path.stem} was annotated '
        f'and saved successfully as {output_path.name}. The gene types annotated:\n'
        '\n'.join(all_gene_types)
    )


if __name__ == '__main__':
    # Get command-line arguments
    args = parse_cl_paths({'--output_path': '-o'})

    # Annotate the .h5 file, saving as "output_path"
    message = annotate(args.output_path)

    print(message)
