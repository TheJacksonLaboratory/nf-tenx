#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { construct_library_csv_content } from './functions.nf'


process run_cellranger_arc_count {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool}/*", mode: "copy"
    memory 1.GB
    cpus 1
    tag "$record.output_id"

    container "library://singlecell/${record.tool}:${record.tool_version}"

    input:
      val record
    output:
      tuple val(record), path("${record.tool}/*"), emit: cellranger_arc_outputs
      tuple val(record), path("${record.tool}/*"), emit: hash_dir

    script:
    localmem = Math.round(task.memory.toGiga() * 0.95)
    lib_csv_content = construct_library_csv_content(record)
    """
    lib_csv_file=${record.output_id}.csv
    cat <<-EOF > \$lib_csv_file
$lib_csv_content
EOF

    cellranger-arc count \
      --id=$record.output_id \
      --reference=$record.reference_path \
      --libraries=\$lib_csv_file \
      --localcores=$task.cpus \
      --localmem=$localmem

    # do rearrange here
    mkdir -p ${record.tool}/_files
    mv ${record.output_id}/_* ${record.tool}/_files
    mv ${record.output_id}/*.tgz ${record.tool}/
    find ${record.output_id}/SC_ATAC_GEX_COUNTER_CS -type f -name "summary.json" -exec mv {} ${record.tool}/summary.json \\;
    mv ${record.output_id}/outs/* ${record.tool}/
    mv \$lib_csv_file ${record.tool}/
    """
}
