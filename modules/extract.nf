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
    tuple val(record), path('*')
    val genes_genome_path

    output:
    tuple val(record), path('*', includeInputs: true), emit: cellranger
    path(genes_genome_path)

    script:
    """
    extract_files.py --ref_genome_dir=${record.reference_path} --annots_dir=${params.gene_annotations_dir} --genes_genome_csv=${genes_genome_path}
    """
}