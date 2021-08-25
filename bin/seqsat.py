#!/usr/bin/env python3
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import argparse

import h5py

import pandas as pd
import numpy as np


def decompress_seq(x, length=16, bits=64):
    nucleotides = np.array(list(b"ACGT"))
    x = np.uint64(x)
    assert length <= (bits / 2 - 1)
    if x & (1 << (bits - 1)):
        return "N" * length
    result = bytearray(length)
    for i in range(length):
        result[(length - 1) - i] = nucleotides[x & np.uint64(0b11)]
        x = x >> np.uint64(2)
    return result.decode("ascii")


def compute_legacy(mole_file, filt_mat_file):
    with h5py.File(mole_file) as f:
        print(list(f))
        genome = f["/genome_ids"][0].decode("ascii")

        gem_group = np.array(f["/gem_group"]).astype(str)
        barcodes = np.array(f["/barcode"]).astype(str)
        reads = np.array(f["/reads"])  # .astype(np.int64)
        umis = reads.astype(bool).astype(np.int8)

        gems = np.char.add("-", gem_group)
        bcs = np.char.add(np.array([decompress_seq(bc) for bc in barcodes]), gems)
        uniq_bcs, uniq_coords, inverse_coords = np.unique(
            bcs, return_index=True, return_inverse=True
        )
        rcs = np.bincount(inverse_coords, reads)
        ucs = np.bincount(inverse_coords, umis)

        seqsat = pd.DataFrame({"umis": ucs, "reads": rcs}, index=bcs[uniq_coords])
        seqsat["saturation"] = (seqsat.reads - seqsat.umis) * 100 / seqsat.reads

    with h5py.File(filt_mat_file) as fin:
        filt_bcs = pd.Index(fin[f"/{genome}/barcodes"]).astype(str)

    return seqsat.reindex(filt_bcs)


def compute_modern(mole_file, filt_mat_file, *args):
    with h5py.File(mole_file) as f:

        gem_group = np.array(f["/gem_group"])
        barcode_idx = np.array(f["/barcode_idx"]).astype(np.int64)
        barcodes = np.array(f["/barcodes"]).astype(str)
        reads = np.array(f["/count"])  # .astype(np.int64)
        umis = reads.astype(bool).astype(np.int8)

        uniq_idxs, uniq_coords = np.unique(barcode_idx, return_index=True)
        gems = np.char.add("-", gem_group[uniq_coords].astype(str))
        bcs = np.char.add(barcodes[uniq_idxs], gems)

        rcs = np.bincount(barcode_idx, reads)[uniq_idxs]
        ucs = np.bincount(barcode_idx, umis)[uniq_idxs]

        seqsat = pd.DataFrame({"umis": ucs, "reads": rcs}, index=bcs)
        seqsat["saturation"] = (seqsat.reads - seqsat.umis) * 100 / seqsat.reads

    with h5py.File(filt_mat_file) as fin:
        filt_bcs = pd.Index(fin["/matrix/barcodes"]).astype(str)

    return seqsat.reindex(filt_bcs)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("molecule_info_path", help="Path to 'molecule_info.h5' file")
    parser.add_argument(
        "filtered_matrix_path", help="Path to 'filtered_*_bc_matrix*.h5' file"
    )
    parser.add_argument("software", help="10X software version")
    parser.add_argument("software_version", help="10X software version")
    parser.add_argument(
        "--outfile",
        default="sequencing_saturation.csv",
        help="Output csv file for saturation metrics.",
    )

    return parser.parse_args()


def main():
    args = parse_args()
    if (args.software == "cellranger") and (args.software_version[0] in "12"):
        seqsat = compute_legacy(
            args.molecule_info_path, args.filtered_matrix_path
        )
    else:
        seqsat = compute_modern(args.molecule_info_path, args.filtered_matrix_path)
    seqsat.to_csv(args.outfile)


if __name__ == "__main__":
    main()
