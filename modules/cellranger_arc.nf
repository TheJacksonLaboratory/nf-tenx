#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { construct_library_csv_content; join_map_items } from './functions.nf'
 

def construct_arc_cli_options(record) {
    options = [:]
    options["--id"] = record.output_id
    options["--reference"] = record.reference_path
    options["--libraries"] = "${record.output_id}.csv"
    options["--description"] = record.sample_name
    return(join_map_items(options))
}


process CELLRANGER_ARC_COUNT {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool_pubdir}/*", mode: "copy"
    tag "$record.output_id"

    // cpus determined by profile
    // memory determined by profile
    // n_reads here counts both ATAC + GEX, so scale down the factor here.
    time { (record.n_reads / 300000000).round(2) * 6.hour * params.time_scale_factor }

    container "library://singlecell/${record.tool}:${record.tool_version}"

    input:
      val record
    output:
      tuple val(record), path("${record.tool_pubdir}/*"), emit: cellranger_arc_outputs
      tuple val(record), path("${record.tool_pubdir}/*"), emit: hash_dir
      tuple val(record), path("${record.tool_pubdir}/*{metrics,summary}*", glob: true), emit: metrics

    script:
    main_options = construct_arc_cli_options(record)
    lib_csv_content = construct_library_csv_content(record)
    localmem = Math.round(task.memory.toGiga() * 0.95)
    """
    lib_csv_file=${record.output_id}.csv
    echo -e '$lib_csv_content' > \$lib_csv_file

    cellranger-arc count $main_options --localcores=$task.cpus --localmem=$localmem

    # do rearrange here
    mkdir -p ${record.tool_pubdir}/_files
    mv \$lib_csv_file ${record.tool_pubdir}/
    mv ${record.output_id}/_* ${record.tool_pubdir}/_files
    mv ${record.output_id}/*.tgz ${record.tool_pubdir}/
    mv ${record.output_id}/outs/* ${record.tool_pubdir}/
    find ${record.output_id}/SC_ATAC_GEX_COUNTER_CS -type f -name "summary.json" -exec mv {} ${record.tool_pubdir}/summary_atac.json \\;
    find ${record.output_id}/SC_ATAC_GEX_COUNTER_CS -type f -name "metrics_summary_json.json" -exec mv {} ${record.tool_pubdir}/summary_gex.json \\;
    """
}
