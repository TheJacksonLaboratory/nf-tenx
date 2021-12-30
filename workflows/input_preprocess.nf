#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/
import org.yaml.snakeyaml.Yaml

include { CHECK_SAMPLESHEET; COUNT_READS } from '../modules/preflight'
include { read_n_reads; collect_fastqs; parse_gt_ids; parse_yaml } from '../modules/samplesheet'
include { construct_output_id; construct_tool_pubdir; } from '../modules/samplesheet'

workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    CHECK_SAMPLESHEET(samplesheet).yml \
        | map { parse_yaml(it) } \
        | flatten() \
        | map { construct_output_id(it) } \
        | map { construct_tool_pubdir(it) } \
        | map { collect_fastqs(it) } \
        | map { parse_gt_ids(it) } \
        | COUNT_READS \
        | map { rec, f -> read_n_reads(rec, f) } \
        | set { main_records }

    emit:
    main_records
    versions = CHECK_SAMPLESHEET.out.versions
}

