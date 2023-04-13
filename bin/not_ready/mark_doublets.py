#!/usr/bin/env python3

import scanpy as sc
import scrublet as scr
from argparse import ArgumentParser
from pathlib import Path

# Initialize parser, add argument, and parse
parser = ArgumentParser()
parser.add_argument("--anndata_path", "-a", type=Path, required=True)
args = parser.parse_args()

# Read anndata object into memory
adata = sc.read_h5ad(args.anndata_path)

# Calculate predicted doublets and doublet scores
scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata.obs["doublet_score"] = doublet_scores
adata.obs["doublet_predicted"] = predicted_doublets

output_path = Path("final_anndata.h5ad")
adata.write_h5ad(output_path)
