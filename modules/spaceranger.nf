#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

def construct_vis_cli_options(record) {
    options = [:]
    options["--id"] = record.output_id
    options["--transcriptome"] = record.reference_path
    options["--sample"] = record.prefixes.join(",")
    options["--fastqs"] = record.fastq_paths.join(",")

    options["image"] = record.image
    options["slide"] = record.slide
    options["area"] = record.area
    options["--description"] = record.sample_name

    if (record.dark_image) { options["--darkimage"] = record.dark_image }
    if (record.manual_alignment) { options["--loupe-alignment"] = record.manual_alignment }
    if (record.no_bam) { options["--no-bam"] = "" }

    return(join_map_items(options))
}


process SPACERANGER_COUNT {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool_pubdir}/*", mode: "copy"
    tag "$record.output_id"
    label "10x_genomics_count"

    container "library://singlecell/${record.tool}:${record.tool_version}"

    input:
      val record
    output:
      tuple val(record), path("${record.tool_pubdir}/*"), emit: spaceranger_outputs
      tuple val(record), path("${record.tool_pubdir}/*"), emit: hash_dir
      tuple val(record), path("${record.tool_pubdir}/*{metrics,summary}*", glob: true), emit: metrics

    script:
    main_options = construct_vis_cli_options(record)
    localmem = Math.round(task.memory.toGiga() * 0.95)
    """
    spaceranger count $main_options --localcores=$task.cpus --localmem=$localmem

    mkdir -p ${record.tool_pubdir}/_files
    mv ${record.output_id}/_* ${record.tool_pubdir}/_files
    mv ${record.output_id}/*.tgz ${record.tool_pubdir}/
    mv ${record.output_id}/outs/* ${record.tool_pubdir}/
    find ${record.output_id}/SPATIAL_RNA_COUNTER_CS/ -type f -name "metrics_summary_json.json" -exec mv {} ${record.tool_pubdir}/summary.json \\;
    find ${record.output_id}/SPATIAL_RNA_COUNTER_CS/ -type f -name "alignment_metrics.json" -exec mv {} ${record.tool_pubdir}/spatial/alignment_summary.json \\;
    """
}
