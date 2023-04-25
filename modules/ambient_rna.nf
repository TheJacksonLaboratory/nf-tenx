#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process FILTER_AMBIENT_RNA {
    tag "$record.output_id"
    publishDir params.pubdir, pattern: "*soupx*", mode: "copy"
    label 'alternative_count'

    input:
    tuple val(record), path('*')

    output:
    tuple val(record), path('*')

    script:
    """
    mkdir outs
    filter_ambient_rna.r
    """
}