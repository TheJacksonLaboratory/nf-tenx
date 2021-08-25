#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process SEQUENCING_SATURATION {
    publishDir "${params.pubdir}/${record.output_id}/${record.tool}", pattern: "*", mode: "copy"
    memory 1.GB
    cpus 1
    tag "$record.output_id"

    container "library://singlecell/python:3.8"

    input:
      tuple val(record), path('*')
    output:
      path 'sequencing_saturation.csv'

    script: // packaged under projectDir/bin
    """
    python seqsat.py molecule*.h5 filtered_*.h5 ${record.tool} ${record.tool_version}
    """
