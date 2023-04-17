#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process VELOCYTO {
    tag "$record.output_id"
    publishDir params.pubdir, pattern: '**/*.loom', mode: 'copy'
    label "velocyto"

    input:
    tuple val(record), path('*'), path('*genes*.gtf', stageAs: 'velocyto_genes.gtf')

    output:
    tuple val(record), path('**/*.loom')

    script:
    """
    ln -s . outs
    velocyto run10x --samtools-threads=${task.cpus} --logic=ObservedSpanning10X . velocyto_genes.gtf
    """
}