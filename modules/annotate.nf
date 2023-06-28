#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process ANNOTATE_MTX {
    tag "$record.output_id-$tool"
    executor 'local'
    publishDir "${params.pubdir}/${record.output_id}/annotations/${tool}", mode: "copy", pattern: "*.h5ad"
    container 'library://singlecell/python3.10-scrna-annotation:latest'

    input:
    tuple val(record), path('*'), val(tool)
  
    output:
    tuple val(record), path('*'), val(tool)

    script:
    """
    annotate.py --record_id=${record.output_id} --tool=${tool} --annots_dir=${params.annots_dir}
    """
}
