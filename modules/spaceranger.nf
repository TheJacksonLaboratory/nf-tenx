#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { join_map_items } from './functions.nf'


params.probe_dir = "${workflow.projectDir}/assets/probe_sets"


def construct_vis_cli_options(record) {
    options = [:]
    options["--id"] = record.output_id
    options["--sample"] = record.prefixes.join(",")
    options["--fastqs"] = record.fastq_paths.join(",")
    options["--transcriptome"] = record.reference_path

    options["--image"] = record.image
    options["--slide"] = record.slide
    options["--area"] = record.area
    options["--description"] = record.sample_name

    if (!record.get("probe_set", "").isEmpty()) {
        options["--probe-set"] = file("${params.probe_dir}/${record.probe_set}", checkIfExists: true)
    }

    if (record.dark_image) { options["--darkimage"] = record.dark_image }
    if (record.color_image) { options["--colorizedimage"] = record.color_image }
    if (record.manual_alignment) { options["--loupe-alignment"] = record.manual_alignment }
    if (record.slide_file) { options["--slidefile"] = record.slide_file }
    if (record.requires_rotation) { options["--reorient-images"] = null }
    if (record.no_bam) { options["--no-bam"] = null }

    return(join_map_items(options))
}


process SPACERANGER_COUNT {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool_pubdir}/*", mode: "copy"
    tag "$record.output_id"
    label "tenx_genomics_count"

    // cpus determined by profile
    // memory determined by profile
    time { (record.n_reads / 300000000).round(2) * 8.hour * params.time_scale_factor }

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

process IMAGE_PROCESS {
    publishDir "${params.pubdir}/${record.output_id}/img", pattern: "*", mode: "copy"
    tag "$record.output_id"
    label "tenx_genomics_count"

    container "library://singlecell/${record.tool}:${record.tool_version}"

    input:
      val(record)
    output:
      path("*")
      tuple val(record), path("*.md5"), emit: img_hashes

    script:
    record.roi_json = record.roi_json ?: ""
    """
    cp ${record.image} .
    if [ ! -z "${record.roi_json}" ]; then
        cp ${record.roi_json} .
    fi
    md5sum * > hashes_img.md5
    """
}
