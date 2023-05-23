#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process VELOCYTO {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/annotations/${tool}", pattern: "*velocyto*.loom", mode: "copy"
    label "velocyto"
    container 'docker://mparikhbroad/velocyto'

    input:
    tuple val(record), path('*', stageAs: 'outs/*'), val(tool)

    output:
    tuple val(record), path('*', includeInputs: true), val(tool)

    script:
    genes_gtf = 'velocyto_genes.gtf'
    """
    velocyto run10x --samtools-threads=${task.cpus} --logic=ObservedSpanning10X . outs/${genes_gtf}
    mv **/*.loom ${record.output_id}_${tool}_velocyto.loom
    rmdir velocyto
    """
}

process TEST {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/annotations/${tool}", pattern: "*", mode: "copy"
    container '/sc/service/analysis/tmp/pipeline_development/nextflow-dev/containers/py_w_loompy.sif'

    input:
    tuple val(record), path('*'), val(tool)

    output:
    tuple val(record), path('*', includeInputs: true), val(tool)

    script:
    """
    test.py
    """
}