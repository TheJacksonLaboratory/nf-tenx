#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process run_cellranger_arc_count {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool}/*", mode: "copy"
    memory 1.GB
    cpus 1
    tag "$record.output_id"

    container "library://singlecell/${record.tool}:${record.tool_version}"

    input:
      val record
    output:
      tuple val(record), path("${record.tool}/*"), emit: cellranger_arc_outputs
      tuple val(record), path("${record.tool}/*"), emit: hash_dir

    script:
    samples = record.prefixes.join(",")
    cells = record.n_cells ?: 6000
    localmem = Math.round(task.memory.toGiga() * 0.95)
    """
    cellranger-arc count \
      --id=$record.output_id \
      --sample=$samples \
      --fastqs=$record.fastq_path \
      --transcriptome=$record.reference_path \
      --expect-cells=$cells \
      --localcores=$task.cpus \
      --localmem=$localmem

    # do rearrange here
    mkdir -p ${record.tool}/_files
    mv ${record.output_id}/_* ${record.tool}/_files
    mv ${record.output_id}/*.tgz ${record.tool}/
    mv ${record.output_id}/outs/* ${record.tool}/
    """
}
