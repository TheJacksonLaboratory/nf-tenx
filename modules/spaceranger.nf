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


process run_spaceranger_count {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool}/*", mode: "copy"
    tag "$record.output_id"

    container "library://singlecell/${record.tool}:${record.tool_version}"

    input:
      val record
    output:
      tuple val(record), path("${record.tool}/*"), emit: spaceranger_outputs
      tuple val(record), path("${record.tool}/*"), emit: hash_dir

    script:
    main_options = construct_vis_cli_options(record)
    localmem = Math.round(task.memory.toGiga() * 0.95)
    """
    spaceranger count $main_options --localcores=$task.cpus --localmem=$localmem

    mkdir -p ${record.tool}/_files
    mv ${record.output_id}/_* ${record.tool}/_files
    mv ${record.output_id}/*.tgz ${record.tool}/
    mv ${record.output_id}/outs/* ${record.tool}/
    find ${record.output_id}/SPATIAL_RNA_COUNTER_CS/ -type f -name "metrics_summary_json.json" -exec mv {} ${record.tool}/summary.json \\;
    find ${record.output_id}/SPATIAL_RNA_COUNTER_CS/ -type f -name "alignment_metrics.json" -exec mv {} ${record.tool}/spatial/alignment_summary.json \\;
    """
}
