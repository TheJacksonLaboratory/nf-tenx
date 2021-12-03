#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process CELLRANGER_MULTI {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool_pubdir}/*", mode: "copy"
    tag "$record.output_id"
    label "10x_genomics_other"

    container "library://singlecell/${record.tool}:${record.tool_version}"

    input:
      val record
    output:
      tuple val(record), path("${record.tool_pubdir}/*"), emit: cellranger_outputs
      tuple val(record), path("${record.tool_pubdir}/*"), emit: hash_dir

    script:
    samples = record.prefixes.join(",")
    cells = record.n_cells ?: 6000
    localmem = Math.round(task.memory.toGiga() * 0.95)
    """
    """
}
