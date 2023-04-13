#!/usr/bin/env python3

import scanpy as sc
from argparse import ArgumentParser
from pathlib import Path

# Create a list of expected command-line arguments and their settings
# for the argument parser
expected_args = [
    {
        "long_name": "--velocyto_loom_path",
        "short_name": "-v",
        "required": True,
        "type": Path,
    },
    {
        "long_name": "--mat_path",
        "short_name": "-m",
        "required": True,
        "type": Path,
    },
]

# Initialize parser and iterate over each expected argument, adding to
# parser, and then parse the arguments
parser = ArgumentParser()

for arg in expected_args:
    parser.add_argument(arg["long_name"], arg["short_name"], type=arg["type"], required=arg["required"])
args = parser.parse_args()

adata = sc.read_10x_h5(args.mat_path)
velo_adata = sc.read_loom(args.velocyto_loom_path)

adata.layers["spliced"] = velo_adata.layers["spliced"].copy()
adata.layers["unspliced"] = velo_adata.layers["unspliced"].copy()
adata.layers["ambiguously_spliced"] = velo_adata.layers["ambiguous"].copy()

output_path = Path("combined_anndata.h5ad")
adata.write_h5ad(output_path)
