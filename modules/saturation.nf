#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process SEQUENCING_SATURATION {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/${record.tool}", pattern: "sequencing_saturation.csv", mode: "copy"
    label "high_mem"

    time { (record.n_reads / 600000000).round(2) * 1.hour * params.time_scale_factor }

    container "library://singlecell/python:3.8"

    input:
      tuple val(record), path('*')
    output:
      file 'sequencing_saturation.csv'
      tuple val(record), path('*'), emit: hash_data

    script: // packaged under projectDir/bin
    """
    seqsat.py molecule*.h5 filtered_*.h5 ${record.tool} ${record.tool_version}
    """
}
