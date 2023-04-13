#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
from argparse import ArgumentParser
from pathlib import Path
import scrublet as scr

# Create a list of expected command-line arguments and their settings
# for the argument parser
expected_args = [
    {
        "long_name": "--data_path",
        "short_name": "-d",
        "required": True,
        "type": Path,
        "nargs": None,
    },
    {
        "long_name": "--ref_df_paths",
        "short_name": "-r",
        "required": True,
        "type": Path,
        "nargs": "+",
    }
]

# Initialize parser and iterate over each expected argument, adding to
# parser, and then parse the arguments
parser = ArgumentParser()

for arg in expected_args:
    parser.add_argument(
        arg["long_name"],
        arg["short_name"],
        type=arg["type"],
        nargs=arg["nargs"],
        required=arg["required"],
    )
args = parser.parse_args()

# Read data into anndata object. The function sc.read() (which should 
# be file-type-agnostic) does not work, so here are some ugly
# if-statements
ext = args.data_path.suffix
if ext == ".h5ad":
    adata = sc.read_h5ad(args.data_path)
elif ext == ".h5":
    adata = sc.read_10x_h5(args.data_path)
elif ext == ".loom":
    adata = sc.read_loom(args.data_path)
else:
    raise ValueError(f"{ext} files are not supported.")

# Iterate over each reference dataframe and each gene type and annotate
for df_path in args.ref_df_paths:
    df = pd.read_pickle(df_path)
    gene_types = [col for col in df.columns if df[col].dtype == bool]
    for gene_type in gene_types:
        adata.var[gene_type] = adata.var["gene_ids"].str.upper().isin(df.loc[df[gene_type], "ensembl_gene_id"])

# Make genes/cells names unique
adata.var_names_make_unique()
adata.obs_names_make_unique()

# Calculate QC metrics. Generate a list of annotated gene types first
annotated_gene_types = [col for col in adata.var.columns if adata.var[col].dtype == bool]
sc.pp.calculate_qc_metrics(adata=adata, qc_vars=annotated_gene_types, inplace=True)

# Calculate predicted doublets and doublet scores
scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata.obs["doublet_score"] = doublet_scores
adata.obs["doublet_predicted"] = predicted_doublets

# Save to anndata object to disk
output_path = Path("final_anndata.h5ad")
adata.write_h5ad(output_path)