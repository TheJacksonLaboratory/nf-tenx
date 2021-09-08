#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { construct_library_csv_content; join_map_items } from './functions.nf'


def construct_atac_cli_options(record) {
    options = [:]
    options["--id"] = record.output_id
    options["--reference"] = record.reference_path
    options["--sample"] = record.prefixes.join(",")
    options["--fastqs"] = record.fastq_paths.join(",")
    options["--description"] = record.sample_name

    return(join_map_items(options))
}


process run_cellranger_atac_count {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool}/*", mode: "copy"
    tag "$record.output_id"

    container "library://singlecell/${record.tool}:${record.tool_version}"

    input:
      val record
    output:
      tuple val(record), path("${record.tool}/*"), emit: cellranger_atac_outputs
      tuple val(record), path("${record.tool}/*"), emit: hash_dir

    script:
    main_options = construct_atac_cli_options(record)
    localmem = Math.round(task.memory.toGiga() * 0.95)
    """
    cellranger-atac count $main_options --localcores=$task.cpus --localmem=$localmem

    # do rearrange here
    mkdir -p ${record.tool}/_files
    mv ${record.output_id}/_* ${record.tool}/_files
    mv ${record.output_id}/*.tgz ${record.tool}/
    mv ${record.output_id}/outs/* ${record.tool}/
    """
}
