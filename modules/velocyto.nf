#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process VELOCYTO {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/annotations/${tool}", pattern: "*.loom", mode: "copy"
    label "velocyto"

    input:
    tuple val(record), path('*', stageAs: 'outs/*'), val(tool)
    path genes_gtf

    output:
    tuple val(record), path('*', includeInputs: true), val(tool)

    script:
    """
    velocyto run10x --samtools-threads=${task.cpus} --logic=ObservedSpanning10X . ${genes_gtf}
    mv **/*.loom ${record.output_id}_velocyto_spliced_mtx_${tool}.loom
    rmdir velocyto
    """
}