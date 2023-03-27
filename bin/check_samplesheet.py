#!/usr/bin/env python

import os
import sys
import logging
import argparse
from pathlib import Path
from string import ascii_letters

import yaml
from yaml import load, dump

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


logging.basicConfig(format="{levelname}: {message}", level=logging.INFO, style="{")


ASSET_DIR = Path(__file__).parent.absolute() / ".." / "assets"


def parse_args():
    desc = "Check samplesheet for missing data"
    epilog = "Example usage:  python check_samplesheet.py <SAMPLESHEET>"

    parser = argparse.ArgumentParser(description=desc, epilog=epilog)
    parser.add_argument("samplesheet", type=Path, help="Input samplesheet file")
    parser.add_argument("outsheet", type=Path, help="Output samplesheet file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Extra logging for use outside pipeline")
    return parser.parse_args()


def print_error(record_id, error, context="Line", context_str=""):
    error_str = f"Samplesheet record {record_id} -> {error}"
    if context != "" and context_str != "":
        error_str = f"Samplesheet record {record_id} -> {error}\n{context.strip()}: '{context_str.strip()}'"
    logging.error(error_str)
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

        self.allowed_fields = {
            "use_undetermined", "lanes",
            "n_cells", 
            "is_nuclei", 
            "design",
            "probe_set", 
            "tags", 
            "no_bam", 
        }

    def check_mutually_exclusive_fields(self, record_id, record, keys1, keys2):
        has_keys1 = any([record.get(key, None) for key in keys1])
        has_keys2 = any([record.get(key, None) for key in keys2])
        if has_keys1 and has_keys2:
            print_error(record_id, f"Keys '{','.join(keys1)}' and '{','.join(keys2)}' are mutually exclusive.")

    def check_needed_field(self, record_id, record, key, is_file=False):
        field = record.get(key, None)
        if not field:
            print_error(record_id, f"Specified assay '{'-'.join(record.get('library_types'))}' must have key '{key}'!")
        if is_file and (not Path(field).exists()):
            print_error(record_id, f"Specified file '{key}'='{field}' does not exist!")
        return field


    def check_unneeded_field(self, record_id, record, key, is_file=False):
        field = record.get(key, None)
        if field:
            print_error(record_id, f"Specified assay '{'-'.join(record.get('library_types'))}' cannot have key '{key}'!")
        return field


    def check_records(self, records):
        for k, record in enumerate(records, start=1):
            if record.get("tool", None) != self.tool:
                continue
            elif record.get("command", None) != self.command:
                continue
            self.check(k, record)
            self.additional_checks(k, record)
    

    def check(self, record_id, record):
        # fields, subclasses modify required fields as needed
        missing = self.required_fields - record.keys()
        if missing:
            print_error(
                record_id, 
                f"Record missing required fields for assay {self.assay_type}: {missing}"
            )

        # check additional fields.
        # this union is done here so that we can modify required fields in subclass inits
        unknown_fields = record.keys() - self.required_fields.union(self.allowed_fields)
        if unknown_fields:
            print_error(
                record_id,
                f"Record has unknown fields: {unknown_fields}"
            )

        # check sample name
        # must not have illegal characters, specifically '+', ' ', ','
        # can't be too long
        sample_name = record.get("sample_name", "")
        allowed_characters = set(ascii_letters + "0123456789_-")
        illegal_characters = set(sample_name) - allowed_characters
        if len(illegal_characters) > 0:
            print_error(
                record_id,
                f"Illegal characters in record sample name: {illegal_characters}"
            )

        if len(sample_name) > 48:
            print_error(record_id, f"Record sample name too long: {sample_name}")

        # library types are allowed
        lts = [
            lt in self.allowed_library_types for lt in record.get("library_types", [])
        ]
        if not all(lts):
            print_error(
                record_id,
                f"Record has invalid library types. Expected {self.allowed_library_types}",
            )

        if not (
            len(record.get("libraries", []))
            == len(record.get("library_types", []))
            == len(record.get("fastq_paths", []))
        ):
            print_error(
                record_id,
                f"Inconsistent number of libraries, library_types, and fastq_paths"
            )

        tv = record.get("tool_version") 
        if tv < self.min_tool_version:
            print_error(
                record_id,
                f"Incompatiable tool version '{tv} < {self.min_tool_version}'"
            )

        lanes = record.get("lanes", [])
        include_undetermined = record.get("use_undetermined", False)
        if include_undetermined and not lanes:
            print_error(
                record_id,
                "FASTQ specification includes 'use_undetermined: true' "
                "but does not include a flowcell lane using 'lanes = []'."
                " This behavior is currently unsupported."
            )

        self.check_fastqs_exist(record_id, record)
        self.check_reference_exists(record_id, record)


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
                record_id,
                f"Record specifies tag(s) we can't find in the catalog: '{missing_tags}'"
            )

    def check_fastqs_exist(self, record_id, record):
        n = 0
        for lib in record["libraries"]:
            for fqp in record["fastq_paths"]:
                fqp = Path(fqp)
                if not fqp.exists():
                    print_error(record_id, f"Fastq path {fqp} doesn't exist!")
                n += len(list(fqp.glob(f"{lib}*")))
        if n < (2 * len(record["libraries"])):
            print_error(record_id,
                (
                    f"Cannot find sufficient fastqs for libraries "
                    f"{record['libraries']} under {record['fastq_paths']}"
                )
            )

    def check_reference_exists(self, record_id, record):
        ref_path = record.get("reference_path", None)
        if ref_path is None:
            return
        self.check_needed_field(record_id, record, "reference_path", is_file=True)

    def check_valid_design(self, record_id, record):
        design = self.check_needed_field(record_id, record, "design")

        for i, (bcs, data) in enumerate(design.items()):
            if '.' in bcs:
                print_error(
                    record_id,
                    f"Design barcodes cannot contain '.' (design item: {i})"
                )
            required_keys = set(["description", "name"])
            if required_keys - set(data.keys()) != set():
                print_error(
                    record_id,
                    "Each element in 'design' must specify 'description' and 'name'. "
                    f"Item {i}: does not."
                )

    def check_for_probeset(self, record_id, record, add_msg=""):
        probe_set = record.get("probe_set", None)
        if probe_set is not None:
            # hack so we can test without nextflow
            probe_loc = ASSET_DIR / "probe_sets"
            probe_set = probe_loc / probe_set
            if not probe_set.exists():
                print_error(
                    record_id,
                    f"Cannot find probe_set specified [{probe_set}] under {probe_loc}"
                )
        elif add_msg:
            print_error(add_msg)


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
                    record_id,
                    f"NonGEX library types but no tag list"
                )
            self.check_tags_exist(record_id, record)

        # if a probeset is specified, make sure it exists
        self.check_for_probeset(record_id, record)


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
        # Differentiate between cellplex and flex here:
        types = ["Multiplexing Capture", "LMO", "TotalSeq-A", "TotalSeq-B", "TotalSeq-C"]
        if set(types) & set(record["library_types"]) == set():
            # We must have a probeset then
            extra_error = (
                f"Record {record_id} must either specify a 'probe_set' or "
                f"have one of the following library types: '{types}'"
            )
            self.check_for_probeset(record_id, record, add_msg=extra_error)
        else:
            # n_cells key must be present
            cells = record.get("n_cells", None)
            if cells is None:
                print_error(record_id, f"Must specify 'n_cells'.")

            # is_nuclei key must be present
            nuclei = record.get("is_nuclei", None)
            if nuclei is None:
                print_error(record_id, f"Must specify 'is_nuclei'.")

        # design key must be present
        self.check_valid_design(record_id, record)

        nongex_lib_types = set(self.allowed_library_types) - set(["Gene Expression", "Multiplexing Capture"])
        has_nongex_libs = nongex_lib_types & set(record["library_types"])
        if has_nongex_libs:
            tag_list = record.get("tags", None)
            if not tag_list:
                print_error(
                    record_id,
                    f"NonGEX library types but no tag list or design"
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
                record_id,
                f"CellRanger-ATAC version < 1.2.0 is not recommended.  Use at least 1.2.0"
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
                record_id,
                f"Must specify both library types: {self.allowed_library_types}"
            )


class VISCountChecker(AssayChecker):
    def __init__(self):
        super().__init__(
            "Visium", "spaceranger", "count", ["Spatial Gene Expression", "CytAssist Gene Expression"]
        )
        for key in ("image", ):
            self.required_fields.add(key)
        for key in (
            "slide", "area", 
            "cyta_image", "dark_image", "color_image", "manual_alignment",
            "slide_file", "requires_rotation", "roi_json"
        ):
            self.allowed_fields.add(key)

    def additional_checks(self, record_id, record):
        image_path = self.check_needed_field(record_id, record, "image", is_file=True)
        self.check_needed_field(record_id, record, "slide")
        self.check_needed_field(record_id, record, "area")

        # CytAssist
        cytaimage_path = record.get("cyta_image", None)
        is_cyta = all(["CytAssist" in lib_type for lib_type in record.get("library_types")])
        if is_cyta:
            self.check_needed_field(record_id, record, "cyta_image", is_file=True)
        else:
            self.check_unneeded_field(record_id, record, "cyta_image")

        # FFPE
        self.check_for_probeset(record_id, record)


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
        logging.debug("YAML parsed successfully")
    except yaml.composer.ComposerError as e:
        logging.error("YAML cannot be parsed. See errors below.")
        logging.error(e)

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
    if args.verbose:
        logging.setLevel(logging.DEBUG)

    records = check_samplesheet(args.samplesheet)
    logging.info("Samplesheet passed all checks.")
    with open(args.outsheet, "w") as fout:
        dump(records, fout, Dumper=Dumper)


if __name__ == "__main__":
    sys.exit(main())
