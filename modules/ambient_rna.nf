#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process FILTER_AMBIENT_RNA {
    tag "$record.output_id"
    label 'alternative_count'
    container '/sc/service/analysis/tmp/pipeline_development/nextflow-dev/containers/r-seurat-base-soupx-dropletutils_latest.sif'

    input:
    tuple val(record), path('*'), val(tool)

    output:
    tuple val(record), path('*', includeInputs: true), val('soupx')

    script:
    soupx_dir = 'soupx'
    genome_csv = 'genes_genome.csv'
    """
    mkdir ${soupx_dir}
    filter_ambient_rna.r --args ${genome_csv} ${soupx_dir}
    """
}

process MIMIC_CELLRANGER {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/annotations", pattern: "soupx/*filtered*", mode: "copy"
    container '/sc/service/analysis/tmp/pipeline_development/nextflow-dev/containers/py_w_loompy.sif'

    input:
    tuple val(record), path('*'), val(tool)

    output:
    tuple val(record), path("${tool}/*", includeInputs: true), val(tool)

    script:
    """
    mimic_cellranger.py --soupx_dir=${tool}
    """
}