#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process EXTRACT_FILES {
    tag "$record.output_id"
    executor 'local'
    container '/sc/service/analysis/tmp/pipeline_development/nextflow-dev/containers/py_w_loompy.sif'
    
    input:
    tuple val(record), path('*')

    output:
    tuple val(record), path('*', includeInputs: true), val('cellranger'), emit: cellranger

    script:
    genome_csv_name = "genes_genome.csv"
    genes_gtf_name = "velocyto_genes.gtf"
    """
    extract_files.py --ref_genome_dir=${record.reference_path} --genome_csv_name=${genome_csv_name} --genes_gtf_name=${genes_gtf_name}
    """
}