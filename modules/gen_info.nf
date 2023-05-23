#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process GEN_SUMMARY {
    tag "$record.output_id"
    executor 'local'
    publishDir "${params.pubdir}/${record.output_id}/annotations/${tool}", mode: "copy"
    container '/sc/service/analysis/tmp/pipeline_development/nextflow-dev/containers/py_w_loompy.sif'
    
    input:
    tuple val(record), path('*'), val(tool)
    path summary_dir

    output:
    tuple val(record), path('*')
    
    script:
    """
    gen_plots.py
    gen_summary.py --summary_dir=${summary_dir} --pubdir=${launchDir / params.pubdir / record.output_id}/annotations/
    """
}
