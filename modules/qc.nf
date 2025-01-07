#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


process FASTQC {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/fastq/qc", pattern: "*", mode: "copy"

    time { (record.n_reads / 300000000).round(2) * 4.hour * params.time_scale_factor }
    cpus 16

    def mem_extract = params.max_memory.toString().replaceAll("[^0-9]", "")
    def int_max_mem = mem_extract.toInteger()
    memory "${ Math.min(int_max_mem,(record.n_reads / 50000000).round(0).toInteger() * 4) }GB" // takes the smallest of the 2 between the max memory set in the config params versus the scaled memory calculation based on read counts

    container "library://singlecell/fastqc:0.11.9"

    input:
      val record

    output:
      tuple val(record), path('*'), emit: fastqc_results

    script:
    files = record.fastqs.join(" ")
    """
    fastqc --threads $task.cpus --quiet --outdir . $files
    """
}


process MULTIQC {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/fastq/qc", pattern: "*", mode: "copy"

    time 4.hour

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

