#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process PREPARE_SEURAT {
    tag "$record.output_id-$tool"
    executor 'local'
    container 'singlecell/python3.10-scrna-annotation:latest'
    
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
    tag "$record.output_id-$tool"
    executor 'local'
    publishDir "${params.pubdir}/${record.output_id}/annotations/${tool}", mode: "copy"
    container 'library://singlecell/r4.2-soupx-seurat:latest'
    
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
