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


def _check_minimal_fields(record_id, record, *additional_fields):
    req_fields = {
        "libraries",
        "library_types",
        "fastq_paths",
        "tool",
        "tool_version",
        "command",
        "sample_name",
        #"reference_path",
    } | set(additional_fields)

    missing = req_fields - record.keys()
    if missing:
        print_error(f"Record {record_id}: missing required fields: {missing}")


def _check_library_types(record_id, record, lib_types):
    lts = [lt in lib_types for lt in record["library_types"]]
    if not all(lts):
        print_error(
            f"Record {record_id}: invalid library types found. Expected {lib_types}",
        )


def _check_n_libs(record_id, record):
    nlibs = len(record["libraries"])
    nlib_types = len(record["library_types"])
    nfastq_paths = len(record["fastq_paths"])
    if not (nlibs == nlib_types == nfastq_paths):
        print_error(
            f"Record {record_id}: inconsistent number of libraries, library_types, and fastq_paths"
        )


def _check_fastqs_exist(record_id, record):
    n = 0
    for lib in record["libraries"]:
        for fqp in record["fastq_paths"]:
            fqp = Path(fqp)
            if not fqp.exists():
                print_error(f"Record {record_id}: fastq path {fqp} doesn't exist!")
            n += len(list(fqp.glob(f"{lib}*")))
    if n < (2 * len(record["libraries"])):
        print_error(
            (
                "Record {record_id}: cannot find sufficient fastqs for libraries "
                f"{record['libraries']} under {record['fastq_paths']}"
            )
        )


def _check_atac_version(record_id, record):
    if record["tool_version"][:3] in ["1.0", "1.1"]:
        print_error(
            (
                "Record {record_id}: CellRanger-ATAC version < 1.2.0 is not recommended.  Use at least 1.2.0"
            )
        )


def _check_image_exists(record_id, record):
    if not Path(record.image).exists():
        print_error("Record {record_id}: image {record.image} does not exist")


def _check_ffpe_probeset(record_id, record):
    # hack so we can test without nextflow
    bin_loc = Path(__file__).parent.absolute()
    assets_loc1 = bin_loc / ".. " / "assets"
    probe_loc = assets_loc1 / "probe_sets"
    probe_set = probe_loc / record.probe_set
    if not probe_set.exists():
        print_error(
            "Record {record_id}: cannot find probe_set specified [{record.probe_set}] under {probe_loc}"
        )


def check_assay(
    assay_name, assay_tool, software_commands, allowed_library_types, records, **kwargs
):
    custom_funcs = kwargs.get("additional_checks", [])
    custom_fields = kwargs.get("additional_fields", [])
    for k, record in enumerate(records, start=1):
        if (record["tool"] != assay_tool) or (record["command"] not in software_commands):
            continue
        _check_minimal_fields(k, record, *custom_fields)
        _check_library_types(k, record, allowed_library_types)
        _check_n_libs(k, record)
        _check_fastqs_exist(k, record)
        for func in custom_funcs:
            func(k, record)


def check_samplesheet(samplesheet):
    # make sure theres no formatting errors
    # this will also ensure multiple documents are not used
    try:
        with open(samplesheet, "r") as fin:
            records = load(fin, Loader=Loader)
    except yaml.composer.ComposerError as e:
        print_error(e)

    check_assay(
        "GEX",
        "cellranger",
        ["count", "", None],
        ["Gene Expression", "Antibody Capture", "CRISPR Guide Capture"],
        records,
    )

    check_assay(
        "VDJ",
        "cellranger",
        ["vdj"],
        ["Immune Profiling", "TCR", "BCR", "VDJ"],
        records,
    )

    check_assay(
        "ATAC",
        "cellranger-atac",
        ["count", "", None],
        ["Chromatin Accessibility"],
        records,
        additional_checks=[_check_atac_version],
    )

    check_assay(
        "ARC",
        "cellranger-arc",
        ["count", "", None],
        ["Gene Expression", "Chromatin Accessibility"],
        records,
    )

    check_assay(
        "Visium",
        "spaceranger",
        ["count", "", None],
        ["Spatial Gene Expression"],
        records,
        additional_checks=[_check_image_exists, _check_ffpe_probeset],
        additional_fields=["slide", "area", "image", "probe_set"],
    )

    check_assay(
        "CITE-SEQ",
        "citeseq_count",
        ["", None],
        ["TotalSeq-A", "TotalSeq-B", "TotalSeq-C", "CMO", "LMO"],
        records,
        additional_fields=["tags"],
    )


def main():
    args = parse_args()
    check_samplesheet(args.samplesheet)
    print("\tSamplesheet passed all checks.")


if __name__ == "__main__":
    sys.exit(main())
