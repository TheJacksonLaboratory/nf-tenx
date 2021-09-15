#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { construct_library_csv_content; join_map_items } from './functions.nf'


def construct_vdj_cli_options(record) {
    options = [:]
    options["--id"] = record.output_id
    options["--reference"] = record.reference_path
    options["--sample"] = record.prefixes.join(",")
    options["--fastqs"] = record.fastq_paths.join(",")
    options["--description"] = record.sample_name
    return(join_map_items(options))
}

process run_cellranger_vdj {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool_pubdir}/*", mode: "copy"
    tag "$record.output_id"

    // cpus determined by profile
    // memory determined by profile
    time { (record.n_reads / 300000000).round(2) * 4.hour * params.time_scale_factor }

    container "library://singlecell/${record.tool}:${record.tool_version}"
 
    input:
      val(record)
    output:
      tuple val(record), path("${record.tool_pubdir}/*"), emit: cellranger_vdj_outputs
      tuple val(record), path("${record.tool_pubdir}/*"), emit: hash_dir
      tuple val(record), path("${record.tool_pubdir}/*{metrics,summary}*", glob: true), emit: metrics

    script:
    main_options = construct_vdj_cli_options(record)
    localmem = Math.round(task.memory.toGiga() * 0.95)
    """
    #cellranger vdj $main_options --localcores=$task.cpus --localmem=$localmem
    #
    ## Do rearrange here
    #mkdir -p ${record.tool_pubdir}/_files
    #mv ${record.output_id}/_* ${record.tool_pubdir}/_files
    #mv ${record.output_id}/*.tgz ${record.tool_pubdir}/ 
    #mv ${record.output_id}/outs/* ${record.tool_pubdir}/
    #find ${record.output_id}/SC_VDJ_ASSEMBLER_CS/ -type f -name "*_summary_json.json" -exec mv {} ${record.tool_pubdir}/summary.json \\;
    """
}
