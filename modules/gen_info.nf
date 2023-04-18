#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process GEN_PLOTS {
    tag "$record.output_id"
    executor 'local'
    publishDir params.pubdir, mode: 'copy'
    
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
    tag "$record.output_id"
    executor 'local'
    publishDir params.pubdir, mode: 'copy'

    input:
    tuple val(record), path('*')
    path summary_dir

    output:
    tuple val(record), path('*')

    script:
    """
    gen_summary.py --summary_dir=${summary_dir} --pubdir=${launchDir / params.pubdir}
    """
}