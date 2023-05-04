#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process FILTER_AMBIENT_RNA {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/annotations", pattern: "soupx/*filtered*", mode: "copy"
    label 'alternative_count'
    
    input:
    tuple val(record), path('*'), val(tool)
    val genome_csv

    output:
    tuple val(record), path("${soupx_dir}/*", includeInputs: true), val('soupx')

    script:
    soupx_dir = 'soupx'
    """
    mkdir ${soupx_dir}
    filter_ambient_rna.r --args ${genome_csv} ${soupx_dir}
    mimic_cellranger.py --soupx_dir=${soupx_dir}
    """
}