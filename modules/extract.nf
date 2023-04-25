#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process EXTRACT_FILES {
    tag "$record.output_id"
    executor 'local'
    //publishDir params.pubdir, pattern: "*.pickle", mode: "copy"

    input:
    tuple val(record), path('*'), path('*')

    output:
    tuple val(record), path('*')

    script:
    """
    extract_files.py --ref_genome_dir=${record.reference_path} --gene_annotations_dir=${params.gene_annotations_dir}
    """
}