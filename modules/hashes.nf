#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process compute_fastq_hashes {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/fastq", pattern: "*", mode: "copy"

    input:
      val record
    output:
      path("*.fastq.gz")
      path('hashes_fastq.md5')

    script:
    files = record.fastqs.join(" ")
    """
    cp $files .
    md5sum $files > hashes_fastq.md5
    """
}

process compute_processed_hashes {
    tag "$record.output_id"
    publishDir "${params.pubdir}/${record.output_id}/${record.tool}", pattern: "*", mode: "copy"

    input:
      tuple val(record), path(hash_dir)
    output:
      path 'hashes_processed.md5'

    script:
    """
    find -L $hash_dir -type f -exec md5sum {} + > hashes_processed.md5
    """
}

