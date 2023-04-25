#!/usr/bin/env python3

import scanpy as sc
from subprocess import run
from pathlib import Path
from scipy.io import mmwrite

# The only thing staged in by Nextflow should be the anndata file (extension .h5ad).
# Use next() to get it out of the glob generator
adata_path = next(Path.cwd().glob('*.h5ad'))
adata = sc.read_h5ad(adata_path)

# Convert AnnData object to Seurat object and save to disk
# Set conversion directory, where the files used to construct seurat object will be stored
conversion_dir = Path('seurat_conversion')
conversion_dir.mkdir()

# Generate a list of directory names based on RNA-velocity.
# This will work whether or not RNA velocity was run because adata.layers maybe empty
sub_dirs = ['total_counts'] + list(adata.layers)
matrix_dirs = [conversion_dir / Path(d) for d in sub_dirs]
for direc in matrix_dirs:
    direc.mkdir()

    # Save the gene names to barcodes.tsv. By doing headers=False and 
    # columns=[], pandas will write only the indices, which are cell barcodes
    barcodes_path = direc / Path('barcodes.tsv')
    adata.obs.to_csv(barcodes_path, header=False, columns=[], sep='\n')

    # Mimic the features.tsv file outputted by cellranger
    features_path = direc / Path('features.tsv')
    with features_path.open(mode='w') as f:
        for gene in adata.var_names:
            f.write(f'{adata.var.loc[gene, "gene_ids"]}\t{gene}\tGene Expression\n')

    # Write the matrix file. Note that it is a transpose per Seurat's requirements
    matrix_path = direc / Path('matrix')
    if direc.name == 'total_counts':
        mmwrite(matrix_path, adata.X.T)
    else:
        mmwrite(matrix_path, adata.layers[direc.name].T, field='integer')

    # Seurat expects gzipped files, so zip everything in matrix files directory
    run(rf'gzip {direc}/*', shell=True)

# Convert boolean columns to integer for ease of use in R
for column in (col for col in adata.obs.columns if adata.obs[col].dtype == bool):
    adata.obs[column] = adata.obs[column].astype(int)

for column in (col for col in adata.var.columns if adata.var[col].dtype == bool):
    adata.var[column] = adata.var[column].astype(int)

# Save gene/cell annotations as CSVs to add to Seurat object. Note that
# index_label=False for easier compatability with R, per the
# pandas documentation.
obs_path = conversion_dir / Path('obs.csv')
adata.obs.to_csv(obs_path, index_label=False)

var_path = conversion_dir / Path('var.csv')
adata.var.to_csv(var_path, index_label=False)
