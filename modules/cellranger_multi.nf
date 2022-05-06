#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { construct_cellplex_library_csv_content; create_feature_reference; join_map_items } from './functions.nf'


def construct_multi_cli_options(record) {
    options = [:]
    options["--id"] = record.output_id

    options["--csv"] = "${record.output_id}.csv"
    options["--disable-ui"] = null

    multi_content = construct_cellplex_library_csv_content(record)

    // handle feature barcoding here
    // --feature-ref
    // Example:
    //   If using TotalSeq-B Abs, you must specify the library as 'Antibody
    //   Capture' for feature barcoding to kick in
    //   Otherwise, you must list it as a separate library with tool cite-seq 
    // must specify "tags" and "feature_ref_type" in record
    feature_ref_types = ["Custom", "Antibody Capture", "CRISPR Guide Capture"]
    ref_content = ""
    if (
        (!record.libraries.intersect(feature_ref_types).isEmpty()) &&
        (record.get("feature_ref_type", "") in feature_ref_types) &&
        (!record.get("tags", []).isEmpty())
    ) {
        multi_content += "\n[feature]\nreference,${record.output_id}_feature_ref.csv\n"
        ref_content = create_feature_reference(record)
    }

    return [join_map_items(options), multi_content, ref_content]
}


process CELLRANGER_MULTI {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool_pubdir}/*", mode: "copy"
    tag "$record.output_id"
    label "tenx_genomics_count"

    container "library://singlecell/${record.tool}:${record.tool_version}"

    time { (record.n_reads / 300000000).round(2) * 8.hour * params.time_scale_factor }

    input:
      val record
    output:
      tuple val(record), path("${record.tool_pubdir}/*"), emit: cellranger_outputs
      tuple val(record), path("${record.tool_pubdir}/*"), emit: hash_dir
      tuple val(record), path("${record.tool_pubdir}/*{metrics,summary}*", glob: true), emit: metrics

    script:
    def (main_options, multiplex_content, feature_ref_content) = construct_multi_cli_options(record)
    def localmem = Math.round(task.memory.toGiga() * 0.95)
    """
    multi_csv_file=${record.output_id}.csv
    echo -e '$multiplex_content' > \$multi_csv_file

    if [ ! -z "${feature_ref_content}" ]; then
        echo -e '$feature_ref_content' > ${record.output_id}_feature_ref.csv
    fi

    cellranger multi $main_options --localcores=$task.cpus --localmem=$localmem

    # do rearrange here
    mkdir -p ${record.tool_pubdir}/_files
    mv \$multi_csv_file ${record.tool_pubdir}/
    mv ${record.output_id}/_* ${record.tool_pubdir}/_files
    mv ${record.output_id}/*.tgz ${record.tool_pubdir}/
    mv ${record.output_id}/outs/* ${record.tool_pubdir}/
    find ${record.output_id}/SC_MULTI_CS/ -type f -name "*_summary_json.json" -exec mv {} ${record.tool_pubdir}/summary.json \\;
    find ${record.output_id}/SC_MULTI_CS/ -type f -name "tag_contaminant_info.json" -exec mv {} ${record.tool_pubdir}/multi/tag_contaminant_info.json \\;
    find ${record.output_id}/SC_MULTI_CS/ -type f -name "marginal_tag_call_metrics_json.json" -exec mv {} ${record.tool_pubdir}/multi/marginal_tag_calls.json \\;
    find ${record.output_id}/SC_MULTI_CS/ -type f -name "non_tag_assignments.json" -exec mv {} ${record.tool_pubdir}/multi/non_tag_assignments.json \\;
    find ${record.output_id}/SC_MULTI_CS/ -type f -name "tag_umi_thresholds_json.json" -exec mv {} ${record.tool_pubdir}/multi/tag_umi_thresholds.json \\;
    """
}
