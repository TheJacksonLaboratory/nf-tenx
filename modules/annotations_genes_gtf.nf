#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process GET_ANNOTATIONS_AND_GENES_GTF {
    tag "$record.output_id"
    executor 'local'
    publishDir params.pubdir, pattern: "*.pickle", mode: "copy"

    input:
    val record

    output:
    tuple val(record), path('*')
    path '*.pickle', emit: gene_annotations_paths
    path 'found_genomes.txt', emit: found_genomes_path

    script:
    """
    get_annotations_and_genes_gtf.py --ref_genome_dir=${record.reference_path} --gene_annotations_dir=${params.gene_annotations_dir}
    """
    }