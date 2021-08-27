#!/usr/bin/env python

import os
import sys
import argparse
from pathlib import Path

import yaml
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


def parse_args():
    desc = "Check samplesheet for missing data"
    epilog = "Example usage:  python check_samplesheet.py <SAMPLESHEET>"

    parser = argparse.ArgumentParser(description=desc, epilog=epilog)
    parser.add_argument("samplesheet", type=Path, help="Input samplesheet file")
    return parser.parse_args()


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)


def _check_library_types(record_id, record, lib_types):
    if not all([lt in lib_types for lt in record["library_types"]]):
        print_error(f"Invalid library types found in record {record_id}", context="Record")


def _check_n_libs(record_id, record):
    nlibs = len(record["libraries"])
    nlib_types = len(record["library_types"])
    nfastq_paths = len(record["fastq_paths"])
    if not (nlibs == nlib_types == nfastq_paths):
        print_error(f"Record {record_id} has inconsistent number of libraries, library_types, and fastq_paths")


def _check_fastqs_exist(record_id, record):
    n = 0
    for lib in record["libraries"]:
        for fqp in record["fastq_paths"]:
            fqp = Path(fqp)
            if not fqp.exists():
                print_error(f"Fastq path {fqp} of record {record_id} doesn't exist!")
            n += len(list(fqp.glob(f"{lib}*")))
    if n < (2 * len(record["libraries"])):
        print_error((
            "Can't find sufficient fastqs for libraries "
            f"{record['libraries']} under {record['fastq_paths']}"
        ))


def check_assay(assay_name, assay_tool, allowed_library_types, records):
    for k, record in enumerate(records, start=1):
        if record['tool'] != assay_tool: continue
        _check_library_types(k, record, allowed_library_types)
        _check_n_libs(k, record)
        _check_fastqs_exist(k, record)


def check_samplesheet(samplesheet):
    # make sure theres no formatting errors
    # this will also ensure multiple documents are not used
    try:
        with open(samplesheet, "r") as fin:
            records = load(fin, Loader=Loader)
    except yaml.composer.ComposerError as e:
        print_error(e)
        
    check_assay("GEX", "cellranger", ["Gene Expression"], records)
    check_assay("ATAC", "cellranger-atac", ["Chromatin Accessibility"], records)
    check_assay("ARC", "cellranger-arc", ["Gene Expression", "Chromatin Accessibility"], records)
    check_assay("Visium", "spaceranger", ["Spatial Gene Expression"], records)
    check_assay("Visium-FFPE", "spaceranger", ["Spatial Gene Expression"], records)
    check_assay("VDJ", "cellranger", ["VDJ Gene Expression"], records)
    
     
def main():
    args = parse_args()
    check_samplesheet(args.samplesheet)
    print("\tSamplesheet passed all checks.")


if __name__ == "__main__":
    sys.exit(main())
