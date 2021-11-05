#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { CHECK_INPUT; COUNT_READS } from '../modules/preflight'
include { read_n_reads; collect_fastqs; parse_gt_ids; load_sample_sheet } from '../modules/samplesheet'
include { construct_output_id; construct_tool_pubdir; } from '../modules/samplesheet'


workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    CHECK_INPUT(samplesheet)

    ch = Channel.fromList(
        load_sample_sheet(CHECK_INPUT.out.yml)
    )
        .map{ it -> construct_output_id(it) }
        .map{ it -> construct_tool_pubdir(it) }
        .map{ it -> collect_fastqs(it) }
        .map{ it -> parse_gt_ids(it) }
         

    COUNT_READS(ch)
        .map { rec, f -> read_n_reads(rec, f) }
        .set { main_records }


    emit:
    main_records
    versions = CHECK_INPUT.out.versions
}

