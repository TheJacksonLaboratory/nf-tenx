#!/usr/bin/env python3

from glob import glob
from pathlib import Path
from shutil import copy

import pandas as pd
import scanpy as sc
import doubletdetection as dd
from arg_utils import parse_cl


def annotate(record_id: str, tool: str, annots_dir: Path) -> str:
    """
    This function takes a filtered feature matrix file and annotates it,
    saving the result as '{record_id}_{tool}_annotated.h5ad'.

    It annotates genes with type and cells with doublet scores (via
    DoubletDetection), also calculating QC metrics while it's at it.
    """
    # Define glob patterns for each matrix type
    filt = '*filtered*matrix*.h5'
    velo = f'*{record_id}*{tool}*velocyto*.loom'

    mtx_readers = {filt: sc.read_10x_h5, velo: sc.read_loom}

    # Read files into a dict of tuples using dict-comprehension and
    # globbing of patterns. Each key is a pattern and each value is a
    # tuple of anndata objects whose file names match that pattern
    data_mtxs = {
        pattern: tuple(reader(path) for path in glob(f'**/{pattern}', recursive=True))
        for pattern, reader in mtx_readers.items()
    }

    # Check that there's only one matrix of each type
    if any(len(path_list) > 1 for path_list in data_mtxs.values()):
        raise ValueError(
            'Something went wrong. The following is a list of glob patterns and a '
            'list of files that matched each glob pattern. At this point of the '
            f'pipeline, only one file should match each pattern:\n{data_mtxs.items()}'
        )

    # Read annotated data object from filtered mtx
    adata = data_mtxs[filt][0]

    # Add velocyto layers to adata if it's there
    if data_mtxs[velo]:
        velo_adata = data_mtxs[velo][0]
        for l in (layer for layer in velo_adata.layers if layer != 'matrix'):
            adata.layers[l] = velo_adata.layers[l].copy()

    # Get the unique genomes in the .h5 file
    filt_mtx_genomes = adata.var['genome'].drop_duplicates()

    # Use generator comprehension to get the annotation files whose stems
    # (each representing a genome) match the genomes in the series created
    # above
    annot_paths = (Path(p) for p in glob(f'{annots_dir}/**', recursive=True))
    annot_paths = (
        path
        for path in annot_paths
        if filt_mtx_genomes.str.match(
            pat=rf'.*{path.name.split(".")[0]}.*', case=False
        ).any()
        and path.is_file()
    )

    # Iterate and copy each annotation to save as outputs of the pipeline
    copied_annots_dir = Path('ref_annotations')
    copied_annots_dir.mkdir()
    for actual_path in annot_paths:
        copy(src=actual_path, dst=copied_annots_dir / actual_path.name)

    # Initialize a set to store all of the gene types in case the
    # reference dataframes have different columns
    all_gene_types = set()

    # Define a dict of pandas readers
    df_readers = {
        '.parquet': pd.read_parquet,
        '.pickle': pd.read_pickle,
        '.csv': pd.read_csv,
        '.feather': pd.read_feather,
    }

    # Iterate over each ref annotation path, read it in, and annotate
    for path in annot_paths:
        ref_df = df_readers[path.suffix](path)

        # The boolean columns of the reference dataframes are the gene
        # annotations, so get them with set comprehension
        current_gene_types = {
            col for col in ref_df.columns if ref_df[col].dtype == bool
        }

        # Add the gene types of the current iteration to the set
        all_gene_types |= current_gene_types

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

    # Calculate predicted doublets and doublet scores, save in obs
    clf = dd.BoostClassifier()
    labels = clf.fit(adata.X).predict()
    scores = clf.doublet_score()

    adata.obs['doublet_score'] = scores
    adata.obs['doublet_predicted'] = labels

    # Save anndata object to h5ad
    output_path = Path(f'{record_id}_{tool}_annotated.h5ad')
    adata.write(output_path)

    return (
        f'The AnnData object associated with {adata.filename} was annotated '
        f'and saved successfully as {output_path.name}. The gene types annotated:\n'
        '\n'.join(all_gene_types)
    )


if __name__ == '__main__':
    # Get command-line arguments
    args = parse_cl(
        ('--record_id', '-r', str), ('--tool', '-t', str), ('--annots_dir', '-a', Path)
    )

    # Annotate the .h5 file
    message = annotate(
        record_id=args.record_id,
        tool=args.tool,
        annots_dir=args.annots_dir,
    )

    print(message)
