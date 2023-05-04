#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process EXTRACT_FILES {
    tag "$record.output_id"
    executor 'local'
    publishDir "${params.pubdir}/${record.output_id}/reference_gene_annotations", pattern: "*.parquet", mode: "copy"

    input:
    tuple val(record), path('*')

    output:
    tuple val(record), path('*', includeInputs: true), val('cellranger'), emit: cellranger
    val(genome_csv_name), emit: genome_csv
    path(genes_gtf_name), emit: genes_gtf

    script:
    genome_csv_name = "genes_genome.csv"
    genes_gtf_name = "velocyto_genes.gtf"
    """
    extract_files.py --ref_genome_dir=${record.reference_path} --annots_dir=${params.annots_dir} --genome_csv_name=${genome_csv_name} --genes_gtf_name=${genes_gtf_name}
    """
}