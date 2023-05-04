#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process ANNOTATE_MTX {
    tag "$record.output_id"
    executor 'local'
    publishDir "${params.pubdir}/${record.output_id}/annotations/${tool}", mode: "copy"

    input:
    tuple val(record), path('*'), val(tool)
  
    output:
    tuple val(record), path('*'), val(tool)

    script:
    """
    annotate.py --output_path=${record.output_id}_annotated_${tool}_mtx.h5ad
    """
}