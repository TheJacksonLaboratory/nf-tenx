#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { construct_library_csv_content; create_feature_reference; join_map_items } from './functions.nf'

def intronic_options(record, opts) {
    major_version = record.tool_version[0].toInteger()

    if (record.is_nuclei == true)  {
        switch (major_version) {
            case { it >= 7 }:
                opts["--include-introns"] = true
                break;

            case { ( it <= 7 ) && (it > 4) }:
                opts["--include-introns"] = null
                break;

            case { it < 4 }:
                if (!(record.reference_path =~ /pre.*rna/)) {
                    throw new Exception('''
                        Library contains nuclei but software version doesn't support --include-introns
                        and the reference path doesn't appear to be a pre-mrna reference
                    ''')
                }
                break;
        }
    } else {
        if (major_version >= 7) {
            opts["--include-introns"] = false
        }
    }
    return options
}


def construct_gex_cli_options(record) {
    options = [:]
    options["--id"] = record.output_id
    options["--transcriptome"] = record.reference_path

    options["--expect-cells"] = record.n_cells ?: 6000
    options["--description"] = record.sample_name
    options["--disable-ui"] = null

    // no bam
    // --no-bam evaluates to true is present at all
    major_version = record.tool_version[0].toInteger()
    if ((record.no_bam) && (major_version >= 5)) {
        options["--no-bam"] = null
    }

    // need to be super careful here
    // --include-introns flag evaluates to true no matter what
    options = intronic_options(record, options)

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
        options["--feature-ref"] = "${record.output_id}_feature_ref.csv"
        ref_content = create_feature_reference(record)
    }

    // As of cellranger 6, they don't encourage using the libraries option for
    // gene expression only.  Therefore, condition its use on whether
    // ref_content is empty or not
    //
    // Additionally, if 'lanes = []' is specified, we can use that here as well.
    // 'Undetermined*.fastq.gz' should already be included in the `prefixes` here.
    if ( ref_content == "" ) {
        options["--sample"] = record.prefixes.join(",")
        options["--fastqs"] = record.fastq_paths.join(",")
        if ( record.get("lanes") && record.get("use_undetermined") ) {
            options["--lanes"] = record.lanes.join(",")
        }
    } else {
        options["--libraries"] = "${record.output_id}.csv"
    }

    return [join_map_items(options), ref_content]
}


process CELLRANGER_COUNT {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool_pubdir}/*", mode: "copy"
    tag "$record.output_id"
    label "tenx_genomics_count"

    container "library://singlecell/${record.tool}:${record.tool_version}"

    // cpus determined by profile
    // memory determined by profile
    time { (record.n_reads / 300000000).round(2) * 8.hour * params.time_scale_factor }

    input:
      val record
    output:
      tuple val(record), path("${record.tool_pubdir}/*"), emit: cellranger_outputs
      tuple val(record), path("${record.tool_pubdir}/*"), emit: hash_dir
      tuple val(record), path("${record.tool_pubdir}/*{metrics,summary}*", glob: true), emit: metrics

    script:
    def (main_options, feature_ref_content) = construct_gex_cli_options(record)
    def lib_csv_content = construct_library_csv_content(record)
    def localmem = Math.round(task.memory.toGiga() * 0.95)
    """
    lib_csv_file=${record.output_id}.csv
    echo -e '$lib_csv_content' > \$lib_csv_file

    if [ ! -z "${feature_ref_content}" ]; then
        echo -e '$feature_ref_content' > ${record.output_id}_feature_ref.csv
    fi

    cellranger count $main_options --localcores=$task.cpus --localmem=$localmem

    # do rearrange here
    mkdir -p ${record.tool_pubdir}/_files
    mv \$lib_csv_file ${record.tool_pubdir}/
    mv ${record.output_id}/_* ${record.tool_pubdir}/_files
    mv ${record.output_id}/*.tgz ${record.tool_pubdir}/
    find ${record.output_id}/SC_RNA_COUNTER_CS/ -type f -name "*_summary_json.json" -exec mv {} ${record.tool_pubdir}/summary.json \\;
    mv ${record.output_id}/outs/* ${record.tool_pubdir}/
    """
}
