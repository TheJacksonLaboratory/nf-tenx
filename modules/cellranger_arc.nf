#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { construct_library_csv_content ; join_map_items } from './functions.nf'


def construct_arc_cli_options(record) {
  options = [:]
  options["--id"] = record.output_id
  options["--reference"] = record.reference_path
  options["--libraries"] = "${record.output_id}.csv"
  options["--description"] = record.sample_name

  [record.disable_cell_annotation, record.create_bam, record.no_bam].each { ele ->
    if (ele == null) {
      return null
    }

    if (ele !instanceof Boolean) {
      throw new Exception("`disable_cell_annotation`, `create_bam`, and `no_bam` must all be boolean values")
    }
  }


  // Define a separate variable just so we get nice editor completions
  def version_parts = record.tool_version.tokenize(".")
  def major_version = version_parts[0].toInteger()
  def minor_version = version_parts[1].toInteger()

  // the option --disable-cell-annotation is only available for cellranger-arc >= 2.2
  def version_is_ge_2_2 = major_version >= 2 && minor_version >= 2

  // Note that we explicitly check that disable_cell_annotation != null instead of a truthiness check
  if (!version_is_ge_2_2 && record.disable_cell_annotation != null) {
    throw new Exception("the option 'disable_cell_annotation' is only available for cellranger-arc >= 2.2")
  }

  if (record.disable_cell_annotation) {
    options["--disable-cell-annotation"] = ""
  }

  def version_is_ge_2_1 = major_version >= 2 && minor_version >= 1
  if (!version_is_ge_2_1 && record.create_bam != null) {
    throw new Exception("the option 'create_bam' is oly available for cellranger-arc >= 2.1")
  }

  if (record.create_bam) {
    options["--create-bam"] = ""
  }

  if (version_is_ge_2_1 && record.no_bam != null) {
    throw new Exception("the option 'no_bam' is only available for cellranger-arc < 2.1")
  }

  if (record.no_bam) {
    options["--no-bam"] = ""
  }

  return (join_map_items(options))
}


process CELLRANGER_ARC_COUNT {
  publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool_pubdir}/*", mode: "copy"
  tag "${record.output_id}"
  label "tenx_genomics_count"

  // cpus determined by profile
  // memory determined by profile
  // n_reads here counts both ATAC + GEX, so scale down the factor here.
  time { (record.n_reads / 300000000).round(2) * 6.hour * params.time_scale_factor }
  module { "${record.tool}/${record.tool_version}" }

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
    echo -e '${lib_csv_content}' > \$lib_csv_file

    cellranger-arc count ${main_options} --localcores=${task.cpus} --localmem=${localmem}

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
