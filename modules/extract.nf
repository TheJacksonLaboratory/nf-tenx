#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process EXTRACT_FILES {
    tag "$record.output_id"
    executor 'local'
    publishDir params.pubdir, pattern: "*.parquet", mode: "copy"

    input:
    tuple val(record), path('*', stageAs: 'outs/*')
    val genome_csv_name
    val genes_gtf_name

    output:
    tuple val(record), path('*', includeInputs: true), emit: extracted
    path(genes_gtf_name), emit: genes_gtf

    script:
    """
    mkdir soupx
    extract_files.py --ref_genome_dir=${record.reference_path} --annots_dir=${params.annots_dir} --genome_csv_name=soupx/${genome_csv_name} --genes_gtf_name=${genes_gtf_name}
    """
}