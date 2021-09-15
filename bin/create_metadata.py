#!/usr/bin/env python
"""
Nextflow DSL2 tool to collect and generate metadata for a sample processed with this
pipeline.

"""
import re
import csv
import json
import argparse
import datetime
from pathlib import Path

import yaml
from yaml import load, dump

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


METADATA_SCHEMA_VERSION = "1.0.0"


# capture record as JSON
# capture GT ID
# capture hashes
# capture date / runner / machine
# capture metrics


def parse_gt_id_from_fastq(r1_fastq_name):
    fastq_pattern = re.compile("(GT[0-9]{2}[-_][0-9]{5})")
    fastq_match = fastq_pattern.search(r1_fastq_name)
    gt_id = None
    if fastq_match:
        gt_id = fastq_match.groups()[0].replace("_", "-")
    return gt_id
    

def create_pipeline_metadata(payload):
    pipeline = {
        "nextflow_version": "$nextflow.version",
        "nextflow_build": "$nextflow.build",
        "nextflow_build": "$nextflow.build",
        "nextflow_script": "$workflow.scriptName",
        "nextflow_script_hash": "$workflow.scriptId",
        "nextflow_workflow_cmdline": "$workflow.commandLine",
        "nextflow_workflow_mnemonic": "$workflow.runName",
        "nextflow_workflow_workdir": "workflow.workDir",
        "nextflow_workflow_submitter": "$workflow.userName",
        "nextflow_workflow_start": "$workflow.start",
    }
    for key in [
        "sampleid",
        "libraryid",
        "gt_release",
        "assay_type",
        "date",
        "genome",
        "cellranger_version",
        "reference_version",
        "is_nuclei",
        "n_cells",
        "read_depth",
        "fraction_flowcell",
        "end_barcoded",
        "chemistry_version"
    ]:
        pipeline["samplesheet_" + key] = payload.get(key, None)

    return pipeline


def create_metadata(args):
    with open(args.record, "r") as fin:
        record = json.load(fin)

    with open(args.workflow, "r") as fin:
        workflow = json.load(fin)

    metadata = {}
    metadata["metadata_generated_date"] = datetime.datetime.isoformat(
        datetime.datetime.now()
    )
    metadata["metadata_schema_version"] = METADATA_SCHEMA_VERSION

    with open(args.record, "r") as fin:
        metadata["record"] = json.load(fin)

    with open(args.workflow, "r") as fin:
        metadata["workflow"] = json.load(fin)

    checksum_func = lambda line: dict(
        zip(("filename", "checksum"), line.strip().split()[::-1])
    )

    metadata["input_checksums"] = []
    for p in args.input_checksums:
        with open(p, "r") as fin:
            metadata["input_checksums"] += sorted(
                list(map(checksum_func, fin)), key=lambda d: d.keys()
            )

    metadata["output_checksums"] = []
    for p in args.output_checksums:
        with open(p, "r") as fin:
            metadata["output_checksums"] += sorted(
                list(map(checksum_func, fin)), key=lambda d: d.keys()
            )

    metadata["metrics"] = []
    for metrics_path in set(args.metrics):
        if not metrics_path.exists(): continue
        if metrics_path.suffix == ".json":
            print("reading json")
            with open(metrics_path, "r") as fin:
                mets = json.load(fin)
        elif metrics_path.suffix == ".csv":
            print("reading csv")
            with open(metrics_path, "r") as fin:
                mets = next(csv.DictReader(fin, delimiter=",", quotechar="\""))
        elif metrics_path.suffix == ".yaml":
            print("reading yaml")
            try:
                with open(metrics_path, "r") as fin:
                    mets = load(fin, Loader=Loader)
            except yaml.composer.ComposerError as e:
                print_error(e)
        else:
            mets = {}

        mets["file"] = str(metrics_path.parts[-1])
        metadata["metrics"].append(mets)

    with open("pipeline-metadata.json", "w") as fout:
        json.dump(metadata, fout)


def parse_args():
    desc = "Collate metadata for an individual record's processing"
    epilog = "Example usage:  python create_metdata.py --record rec.json --workflow wf.json --out metadata.json"

    parser = argparse.ArgumentParser(description=desc, epilog=epilog)
    parser.add_argument("--record", required=True, type=Path, help="JSON file containing record information")
    parser.add_argument("--workflow", required=True, type=Path, help="JSON file containing workflow information")
    parser.add_argument("--out", type=Path, default="pipeline-metadata.json", help="JSON file containing compiled metadata")
    parser.add_argument("--input-checksums", required=True, type=Path, action="append")
    parser.add_argument("--output-checksums", required=True, type=Path, action="append")
    parser.add_argument("--metrics", required=True, type=Path, nargs="+")

    return parser.parse_args()


def main():
    args = parse_args()
    create_metadata(args)


if __name__ == "__main__":
    main()
