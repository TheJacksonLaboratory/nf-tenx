#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process CONVERT_TO_SEURAT {
    tag "$record.output_id"
    executor 'local'
    publishDir "${params.pubdir}/${record.output_id}/annotations/${tool}", mode: "copy", pattern: "*.rds"
    
    input:
    tuple val(record), path('*'), val(tool)
    
    output:
    tuple val(record), path('*')
    
    script:
    total_counts_dir = 'total_counts'
    """
    prepare_seurat.py --total_counts_dir=${total_counts_dir}
    convert_to_seurat.r --args ${total_counts_dir} ${record.output_id}_annotated_${tool}_mtx.rds
    """
}
