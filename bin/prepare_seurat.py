#!/usr/bin/env python3

import gzip
from io import BytesIO
from pathlib import Path

import scanpy as sc
from arg_utils import parse_cl_paths
from scipy.io import mmwrite


def prepare_seurat(total_counts_dir: Path) -> None:
    # The only thing staged in by Nextflow should be the anndata file (extension .h5ad).
    # Use next() to get it out of the glob generator
    adata_path = next(Path().rglob('*.h5ad'))
    adata = sc.read_h5ad(adata_path)

    # Set a layer 'total_counts' equal to adata.X for easier iteration
    adata.layers[total_counts_dir.name] = adata.X  # type: ignore

    # For every layer of adata, the following three files must exist for seurat to work
    file_types = ('matrix.mtx.gz', 'barcodes.tsv.gz', 'features.tsv.gz')

    # Use double generator comprehension to make paths for each file type in each layer
    files_by_layer = (
        (Path(f'{layer}/{f}') for f in file_types) for layer in adata.layers
    )

    # Iterate over each triplet and set up the directory
    for mtx_path, barcodes_path, features_path in files_by_layer:
        # First, make the directory for this layer
        mtx_path.parent.mkdir()

        # Save the gene names to barcodes.tsv. By doing headers=False and
        # columns=[], pandas will write only the indices, which are
        # cell barcodes
        adata.obs.to_csv(barcodes_path, header=False, columns=[], sep='\n')

        # Mimic the features.tsv file outputted by cellranger
        adata.var[['gene_ids', 'feature_types']].to_csv(
            features_path, sep='\t', header=False
        )

        # Write the matrix file into a bytes buffer for gzip compression
        # Note the transpose, since Seurat expects a transposed matrix
        bytes_buffer = BytesIO()
        mmwrite(bytes_buffer, adata.layers[mtx_path.parent.name].T, field='integer')

        with gzip.open(mtx_path, mode='wb') as f:
            f.write(bytes_buffer.getvalue())

    for df_type in ('obs', 'var'):
        # Get either the obs or var dataframe from the AnnData object
        df = getattr(adata, df_type)

        # Convert boolean columns to integer columns for ease of use in R
        df.loc[:, df.dtypes == bool] = df.loc[:, df.dtypes == bool].astype(int)

        # Save annotations as CSV for Seurat. index_label = False for
        # compatability with R (per pandas documentation)
        df.to_csv(Path(f'{df_type}.csv'), index_label=False)


if __name__ == '__main__':
    args = parse_cl_paths({'--total_counts_dir': '-t'})
    prepare_seurat(args.total_counts_dir)
