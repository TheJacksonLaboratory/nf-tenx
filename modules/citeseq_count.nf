#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { join_map_items } from './functions.nf'


params.tag_list = file("${workflow.projectDir}/assets/tags.csv", checkIfExists: true)


def construct_citeseq_cli_options(record) {
    options = [:]

    options["-R1"] = (record.fastqs.findAll { it =~ /_R1_/ }).sort().join(",")
    options["-R2"] = (record.fastqs.findAll { it =~ /_R2_/ }).sort().join(",")
    options["-o"] = "${record.tool_pubdir}"
    options["-t"] = "${record.output_id}_tags.csv"

    options["-cbf"] = record.cbf ?: 1
    options["-cbl"] = record.cbl ?: 16
    options["-umif"] = record.umif ?: 17
    options["-umil"] = record.umil ?: 28

    options["--bc_collapsing_dist"] = record.bc_collapsing_dist ?: 1
    options["--umi_collapsing_dist"] = record.umi_collapsing_dist ?: 2
    options["--expected_cells"] = record.n_cells ?: 40000
    options["--max-errors"] = record.max_errors ?: 1

    // don't specifiy -T (threads) or offset here

    return(join_map_items(options))
}


def create_tag_reference(record) {
    content = []
    offset = 0
    params.tag_list.eachLine { line ->
        def (tag_type, tag_id, read_num, offset5p, tag_sequence, tag_name) = line.split(",")
        if (tag_type in record["library_types"] && tag_id in record["tags"]) {
            content.add("${tag_sequence},${tag_id}_${tag_name}")
            offset = offset5p
        }
    }
    return [content.join("\n"), offset]
}


process CITESEQ_COUNT {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool_pubdir}/*", mode: "copy"
    tag "$record.output_id"

    container "library://singlecell/${record.tool}:${record.tool_version}"

    // cpus determined by profile
    // memory determined by profile
    time { (record.n_reads / 30000000).round(2) * 8.hour * params.time_scale_factor }

    input:
      val record
    output:
      tuple val(record), path("${record.tool_pubdir}/*"), emit: hash_dir
      tuple val(record), path("${record.tool_pubdir}/run_report.yaml"), emit: metrics

    script:
    def (tag_file_content, trim_offset) = create_tag_reference(record)
    main_options = construct_citeseq_cli_options(record)
    """
    echo -e '$tag_file_content' > ${record.output_id}_tags.csv
    mkdir -p ${record.tool_pubdir}

    CITE-seq-Count $main_options -trim=$trim_offset -T=$task.cpus
    echo -e 'CITE-seq-Count $main_options -trim=$trim_offset -T=$task.cpus' > ${record.tool_pubdir}/run_command

    mv ${record.output_id}_tags.csv ${record.tool_pubdir}/
    """
}
