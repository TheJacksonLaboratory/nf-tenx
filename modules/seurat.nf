#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process PREPARE_SEURAT {
    tag "$record.output_id"
    executor 'local'
    publishDir params.pubdir
    
    input:
    tuple val(record), path('*')
    
    output:
    tuple val(record), path('*')
    
    script:
    """
    prepare_seurat.py
    """
}

process CONVERT_TO_SEURAT {
    tag "$record.output_id"
    executor 'local'
    publishDir params.pubdir

    input:
    tuple val(record), path('*')

    output:
    tuple val(record), path('*')

    script:
    """
    convert_to_seurat.r --args ${params.calc_rna_velo}
    """
}