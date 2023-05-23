#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process PREPARE_SEURAT {
    tag "$record.output_id"
    executor 'local'
    container '/sc/service/analysis/tmp/pipeline_development/nextflow-dev/containers/py_w_loompy.sif'
    
    input:
    tuple val(record), path('*'), val(tool)
    
    output:
    tuple val(record), path('*'), val(tool)
    val total_counts_dir, emit: total_counts_dir
    
    script:
    total_counts_dir = "${record.output_id}_${tool}"
    """
    prepare_seurat.py --total_counts_dir=${total_counts_dir}
    """
}

process CONVERT_TO_SEURAT {
    tag "$record.output_id"
    executor 'local'
    publishDir "${params.pubdir}/${record.output_id}/annotations/${tool}", mode: "copy"
    container '/sc/service/analysis/tmp/pipeline_development/nextflow-dev/containers/r-seurat-base-soupx-dropletutils_latest.sif'
    
    input:
    tuple val(record), path('*'), val(tool)
    val total_counts_dir
    
    output:
    tuple val(record), path('*')
    
    script:
    """
    convert_to_seurat.r --args ${total_counts_dir}
    """
}
