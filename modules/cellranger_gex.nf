#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { construct_library_csv_content } from './functions.nf'


process run_cellranger_count {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool}/*", mode: "copy"
    memory 1.GB
    cpus 1
    tag "$record.output_id"

    container "library://singlecell/${record.tool}:${record.tool_version}"

    input:
      val record
    output:
      tuple val(record), path("${record.tool}/*"), emit: cellranger_outputs
      tuple val(record), path("${record.tool}/*"), emit: hash_dir

    script:
    samples = record.prefixes.join(",")
    fastqs = record.fastq_paths.join(",")
    cells = record.n_cells ?: 6000
    lib_csv_content = construct_library_csv_content(record)
    localmem = Math.round(task.memory.toGiga() * 0.95)
    if ( record.tool_version[0].toInteger() < 4 ) 
        """
        cellranger count \
          --id=$record.output_id \
          --sample=$samples \
          --fastqs=$fastqs \
          --transcriptome=$record.reference_path \
          --expect-cells=$cells \
          --localcores=$task.cpus \
          --localmem=$localmem

        # do rearrange here
        mkdir -p ${record.tool}/_files
        mv ${record.output_id}/_* ${record.tool}/_files
        mv ${record.output_id}/*.tgz ${record.tool}/
        mv ${record.output_id}/SC_RNA_COUNTER_CS/*/*/*/*/*/metrics_summary_json.json ${record.tool}/summary.json
        mv ${record.output_id}/outs/* ${record.tool}/
        """

    else
        """
        lib_csv_file=${record.output_id}.csv
        cat <<-EOF > \$lib_csv_file
$lib_csv_content
EOF
        cellranger count \
          --id=$record.output_id \
          --transcriptome=$record.reference_path \
          --libraries=\$lib_csv_file \
          --expect-cells=$cells \
          --localcores=$task.cpus \
          --localmem=$localmem

        # do rearrange here
        mkdir -p ${record.tool}/_files
        mv ${record.output_id}/_* ${record.tool}/_files
        mv ${record.output_id}/*.tgz ${record.tool}/
        mv ${record.output_id}/SC_RNA_COUNTER_CS/*/*/*/*/*/metrics_summary_json.json ${record.tool}/summary.json
        mv ${record.output_id}/outs/* ${record.tool}/
        mv \$lib_csv_file ${record.tool}/
        """
}
