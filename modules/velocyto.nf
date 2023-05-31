#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process VELOCYTO {
    tag "$record.output_id-$tool"
    publishDir "${params.pubdir}/${record.output_id}/annotations/${tool}", pattern: "*velocyto*.loom", mode: "copy"
    label "velocyto"
    container 'docker://mparikhbroad/velocyto'

    input:
    tuple val(record), path('*', stageAs: 'outs/*'), val(tool)

    output:
    tuple val(record), path('*', includeInputs: true), val(tool)

    script:
    """
    velocyto run10x --samtools-threads=${task.cpus} --logic=ObservedSpanning10X . ${record.reference_path}/**/genes.gtf
    mv **/*.loom ${record.output_id}_${tool}_velocyto.loom
    rmdir velocyto
    """
}
