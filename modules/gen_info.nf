#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process GEN_PLOTS {
    executor 'local'
    publishDir params.pubdir
    
    input:
    tuple val(record), path('*')
    
    output:
    tuple val(record), path('*')
    
    script:
    """
    gen_plots.py
    """
}

process GEN_SUMMARY {
    executor 'local'
    publishDir params.pubdir

    input:
    tuple val(record), path('*'), path('*'), path('*')

    output:
    tuple val(record), path('*')

    script:
    """
    gen_summary.py
    """
}