#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


process COMPUTE_FASTQ_HASHES {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/fastq", pattern: "*", mode: "copy"

    input:
      val(record)
    output:
      path("*.fastq.gz")
      path('hashes_fastq.md5')
      tuple val(record), path("*.md5"), emit: input_hashes

    script:
    files = record.fastqs.join(" ")
    """
    cp $files .
    md5sum $files > hashes_fastq.md5
    """
}

process COMPUTE_PROCESSED_HASHES {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/${record.tool_pubdir}", pattern: "*", mode: "copy"

    input:
      tuple val(record), path(hash_dir)
    output:
      path 'hashes_processed.md5'
      tuple val(record), path("*.md5"), emit: output_hashes

    script:
    """
    find -L $hash_dir -type f -exec md5sum {} + > hashes_processed.md5
    """
}

