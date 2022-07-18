#!/usr/bin/env python

import os
import sys
import argparse
from pathlib import Path
from string import ascii_letters

import yaml
from yaml import load, dump

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


ASSET_DIR = Path(__file__).parent.absolute() / ".." / "assets"


def parse_args():
    desc = "Check samplesheet for missing data"
    epilog = "Example usage:  python check_samplesheet.py <SAMPLESHEET>"

    parser = argparse.ArgumentParser(description=desc, epilog=epilog)
    parser.add_argument("samplesheet", type=Path, help="Input samplesheet file")
    parser.add_argument("outsheet", type=Path, help="Output samplesheet file")
    return parser.parse_args()


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)


class AssayChecker:
    supported_tools = {
        "cellranger": ["count", "vdj", "multi", "aggr"],
        "cellranger-arc": ["count", "aggr"],
        "cellranger-atac": ["count", "aggr"],
        "spaceranger": ["count", "aggr"],
        "citeseq-count": [None],
    }

    def __init__(self, assay_type, tool, command, library_types, min_version='0.0.0'):
        self.assay_type = assay_type

        if tool not in self.supported_tools:
            raise Exception(f"Tool {tool} not yet support by this pipeline")
        self.tool = tool
        self.min_tool_version = min_version

        if command not in self.supported_tools[tool]:
            raise Exception(
                f"{tool} command '{command}' not yet support by this pipeline"
            )
        self.command = command

        self.allowed_library_types = library_types

        self.required_fields = {
            "libraries",
            "library_types",
            "fastq_paths",
            "tool",
            "tool_version",
            "command",
            "sample_name",
            "reference_path",
        }

    def check_records(self, records):
        for k, record in enumerate(records, start=1):
            if record.get("tool", None) != self.tool:
                continue
            elif record.get("command", None) != self.command:
                continue
            self.check(k, record)
            self.additional_checks(k, record)

    def check_fastqs_exist(self, record_id, record):
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

    def check(self, record_id, record):
        # fields, subclasses modify required fields as needed
        missing = self.required_fields - record.keys()
        if missing:
            print_error(
                f"Record {record_id} missing required fields for assay {self.assay_type}: {missing}"
            )

        # check sample name
        # must not have illegal characters, specifically '+', ' ', ','
        # can't be too long
        sample_name = record.get("sample_name", "")
        allowed_characters = set(ascii_letters + "0123456789_-")
        illegal_characters = set(sample_name) - allowed_characters
        if len(illegal_characters) > 0:
            print_error(
                f"Record {record_id} illegal characters in sample name: {illegal_characters}"
            )

        if len(sample_name) > 48:
            print_error(f"Record {record_id} sample name too long: {sample_name}")

        # library types are allowed
        lts = [
            lt in self.allowed_library_types for lt in record.get("library_types", [])
        ]
        if not all(lts):
            print_error(
                f"Record {record_id} has invalid library types. Expected {self.allowed_library_types}",
            )

        if not (
            len(record.get("libraries", []))
            == len(record.get("library_types", []))
            == len(record.get("fastq_paths", []))
        ):
            print_error(
                f"Record {record_id} has inconsistent number of libraries, library_types, and fastq_paths"
            )

        tv = record.get("tool_version") 
        if tv < self.min_tool_version:
            print_error(
                f"Record {record_id} has an incompatiable tool version '{tv} < {self.min_tool_version}'"
            )

        lanes = record.get("lanes", [])
        include_undetermined = record.get("use_undetermined", False)
        if include_undetermined and not lanes:
            print_error(
                f"Record {record_id} FASTQ specification includes"
                "'use_undetermined: true' but does not include a flowcell lane using 'lanes = []'."
                " This behavior is currently unsupported."
            )

        self.check_fastqs_exist(record_id, record)


    def check_tags_exist(self, record_id, record):
        tag_catalog = {}
        tag_catalog_file = ASSET_DIR / "tags.csv"
        with open(tag_catalog_file) as fin:
            fin.readline()
            for line in fin:
                items = line.strip().split(",")
                tag_catalog[items[1]] = items

        missing_tags = set(record["tags"]) - set(tag_catalog.keys())
        if len(missing_tags) > 0:
            print_error(
                f"Record {record_id} has tag(s) we can't find in the catalog: '{missing_tags}'"
            )

    def additional_checks(self, record_id, record):
        pass


class GEXCountChecker(AssayChecker):
    def __init__(self):
        super().__init__(
            "GEX",
            "cellranger",
            "count",
            [
                "Gene Expression",
                "Custom",
                "Antibody Capture",
                "CRISPR Guide Capture",
                "TotalSeq-A",
                "TotalSeq-B",
                "TotalSeq-C",
            ],
        )


    def additional_checks(self, record_id, record):
        nongex_lib_types = set(self.allowed_library_types) - set(["Gene Expression"])
        has_nongex_libs = nongex_lib_types & set(record["library_types"])
        if has_nongex_libs:
            tag_list = record.get("tags", None)
            if not tag_list:
                print_error(
                    f"Record {record_id} has nonGEX library types but no tag list"
                )
            self.check_tags_exist(record_id, record)


class GEXMultiChecker(AssayChecker):
    def __init__(self):
        super().__init__(
            "GEX",
            "cellranger",
            "multi",
            [
                "Gene Expression",
                "Multiplexing Capture",
                "Antibody Capture",
                "CRISPR Guide Capture",
                "LMO",
                "TotalSeq-A",
                "TotalSeq-B",
                "TotalSeq-C",
                "Custom"
            ],
        )


    def additional_checks(self, record_id, record):
        types = ["Multiplexing Capture", "LMO", "TotalSeq-A", "TotalSeq-B", "TotalSeq-C"]
        if set(types) & set(record["library_types"]) == set():
            print_error(
                f"Record {record_id} must have one of the following library types: '{types}'"
            )

        # design key must be present
        design = record.get("design", None)
        if design is None:
            print_error(
                f"Record {record_id} must have a 'design' specification"
            )

        # n_cells key must be present
        cells = record.get("n_cells", None)
        if cells is None:
            print_error(
                f"Record {record_id} must have a 'n_cells' specification"
            )

        nongex_lib_types = set(self.allowed_library_types) - set(["Gene Expression", "Multiplexing Capture"])
        has_nongex_libs = nongex_lib_types & set(record["library_types"])
        if has_nongex_libs:
            tag_list = record.get("tags", None)
            if not tag_list:
                print_error(
                    f"Record {record_id} has nonGEX library types but no tag list or design"
                )
            self.check_tags_exist(record_id, record)
        
        


class ImmuneVDJChecker(AssayChecker):
    def __init__(self):
        super().__init__(
            "VDJ",
            "cellranger",
            "vdj",
            ["Immune Profiling", "TCR", "BCR", "VDJ"],
        )


class ATACCountChecker(AssayChecker):
    def __init__(self):
        super().__init__(
            "ATAC", "cellranger-atac", "count", ["Chromatin Accessibility"]
        )

    def additional_checks(self, record_id, record):
        if record["tool_version"][:3] in ["1.0", "1.1"]:
            print_error(
                f"Record {record_id}: CellRanger-ATAC version < 1.2.0 is not recommended.  Use at least 1.2.0"
            )


class ARCCountChecker(AssayChecker):
    def __init__(self):
        super().__init__(
            "ARC",
            "cellranger-arc",
            "count",
            ["Gene Expression", "Chromatin Accessibility"],
        )

    def additional_checks(self, record_id, record):
        # assert we have both library types
        has_both_lib_types = len(set(record["library_types"])) == 2
        if not has_both_lib_types:
            print_error(
                f"Record {record_id} needs to have both library types: {self.allowed_library_types}"
            )


class VISCountChecker(AssayChecker):
    def __init__(self):
        super().__init__(
            "Visium", "spaceranger", "count", ["Spatial Gene Expression"]
        )
        for key in ("slide", "area", "image"):
            self.required_fields.add(key)

    def additional_checks(self, record_id, record):
        image_path = record.get("image", None)
        if not Path(image_path).exists():
            print_error(f"Record {record_id}: image {record['image']} does not exist")

        # FFPE
        probe_set = record.get("probe_set", None)
        if probe_set is not None:
            # hack so we can test without nextflow
            probe_loc = ASSET_DIR / "probe_sets"
            probe_set = probe_loc / probe_set
            if not probe_set.exists():
                print_error(
                    "Record {record_id}: cannot find probe_set specified [{probe_set}] under {probe_loc}"
                )


class CITESeqChecker(AssayChecker):
    def __init__(self):
        super().__init__(
            "CITEseq",
            "citeseq-count",
            None,
            ["TotalSeq-A", "TotalSeq-B", "TotalSeq-C", "CMO", "LMO"],
        )
        self.required_fields.remove("reference_path")
        self.required_fields.add("tags")

    def additional_checks(self, record_id, record):
        self.check_tags_exist(record_id, record)


def check_samplesheet(samplesheet):
    # make sure theres no formatting errors
    # this will also ensure multiple documents are not used
    try:
        with open(samplesheet, "r") as fin:
            records = load(fin, Loader=Loader)
    except yaml.composer.ComposerError as e:
        print_error(e)

    checkers = [
        GEXCountChecker(),
        GEXMultiChecker(),
        ImmuneVDJChecker(),
        ATACCountChecker(),
        ARCCountChecker(),
        VISCountChecker(),
        CITESeqChecker(),
    ]

    for checker in checkers:
        checker.check_records(records)
    return records


def main():
    args = parse_args()
    records = check_samplesheet(args.samplesheet)
    print("\tSamplesheet passed all checks.")
    with open(args.outsheet, "w") as fout:
        dump(records, fout, Dumper=Dumper)


if __name__ == "__main__":
    sys.exit(main())
