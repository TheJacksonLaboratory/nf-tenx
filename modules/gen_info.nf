#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process GEN_PLOTS {
    tag "$record.output_id"
    executor 'local'
    publishDir "${params.pubdir}/${record.output_id}/annotations/${tool}", mode: "copy"
    
    input:
    tuple val(record), path('*'), val(tool)
    
    output:
    tuple val(record), path('*'), val(tool)
    
    script:
    """
    gen_plots.py
    """
}

process GEN_SUMMARY {
    tag "$record.output_id"
    executor 'local'
    publishDir "${params.pubdir}/${record.output_id}/annotations/${tool}", mode: "copy"

    input:
    tuple val(record), path('*', stageAs: 'plots/*'), val(tool)
    path summary_dir

    output:
    tuple val(record), path('*')

    script:
    plots_dir = 'plots'
    """
    gen_summary.py --summary_dir=${summary_dir} --pubdir=${launchDir / params.pubdir} --plots_dir=${plots_dir}
    """
}