#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process ANNOTATE_WITH_VELO {
    tag "$record.output_id"
    executor 'local'
    publishDir params.pubdir

    input:
    tuple val(record), path('*'), path('*')
  
    output:
    tuple val(record), path('*final*anndata*.h5ad')

    script:
    """
    annotate.py
    """
}

process ANNOTATE_NO_VELO {
    tag "$record.output_id"
    executor 'local'
    publishDir params.pubdir

    input:
    tuple val(record), path('*')
  
    output:
    tuple val(record), path('*final*anndata*.h5ad')

    script:
    """
    annotate.py
    """
}