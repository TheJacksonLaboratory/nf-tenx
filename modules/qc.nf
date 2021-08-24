#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


process FASTQC {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/fastq/qc", pattern: "*", mode: "copy"

    container "library://singlecell/fastqc:0.11.9"

    input:
      val record

    output:
      tuple val(record), path('*'), emit: fastqc_results

    script:
    files = record.fastqs.join(" ")
    """
    fastqc -o . $files
    """
}


process MULTIQC {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/fastq/qc", pattern: "*", mode: "move"

    container "library://singlecell/multiqc:1.11"

    input:
      tuple val(record), path('*')

    output:
      path '*'

    script:
    """
    multiqc .
    """
}

