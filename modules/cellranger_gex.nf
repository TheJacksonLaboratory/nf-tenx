#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { construct_library_csv_content; join_map_items } from './functions.nf'


def construct_gex_cli_options(record) {
    options = [:]
    options["--id"] = record.output_id
    options["--transcriptome"] = record.reference_path

    if ( record.tool_version[0].toInteger() < 4 ) {
        options["--samples"] = record.prefixes.join(",")
        options["--fastqs"] = record.fastq_paths.join(",")
    } else {
        options["--libraries"] = "${record.output_id}.csv"
    }

    options["--expect-cells"] = record.n_cells ?: 6000
    options["--include-introns"] = record.is_nuclei ?: false
    options["--description"] = record.sample_name

    return(join_map_items(options))
}


process run_cellranger_count {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool}/*", mode: "copy"
    tag "$record.output_id"

    container "library://singlecell/${record.tool}:${record.tool_version}"

    input:
      val record
    output:
      tuple val(record), path("${record.tool}/*"), emit: cellranger_outputs
      tuple val(record), path("${record.tool}/*"), emit: hash_dir

    script:
    main_options = construct_gex_cli_options(record)
    lib_csv_content = construct_library_csv_content(record)
    localmem = Math.round(task.memory.toGiga() * 0.95)
    println main_options
    """
    lib_csv_file=${record.output_id}.csv
    cat <<-EOF > \$lib_csv_file
$lib_csv_content
EOF

    cellranger count $main_options --localcores=$task.cpus --localmem=$localmem

    # do rearrange here
    mkdir -p ${record.tool}/_files
    mv ${record.output_id}/_* ${record.tool}/_files
    mv ${record.output_id}/*.tgz ${record.tool}/
    find ${record.output_id}/SC_RNA_COUNTER_CS/ -type f -name "*_summary_json.json" -exec mv {} ${record.tool}/summary.json \\;
    mv ${record.output_id}/outs/* ${record.tool}/
    """
}
