#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { join_map_items } from './functions.nf'



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

    if (record.cyta_image) { options["--cytaimage"] = record.cyta_image }
    if (record.dark_image) { options["--darkimage"] = record.dark_image }
    if (record.color_image) { options["--colorizedimage"] = record.color_image }
    if (record.manual_alignment) { options["--loupe-alignment"] = record.manual_alignment }
    if (record.slide_file) { options["--slidefile"] = record.slide_file }
    if (record.requires_rotation) { options["--reorient-images"] = null }

    major_version = record.tool_version[0].toInteger()
    if ((record.no_bam) && (major_version < 3)) { options["--no-bam"] = null }
    if ((major_version >= 3)) { options["--create-bam"] = record.no_bam ? false : true }

    return(join_map_items(options))
}


process SPACERANGER_COUNT {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "${record.tool_pubdir}/*", mode: "copy"
    tag "$record.output_id"
    label "tenx_genomics_count"

    // cpus determined by profile
    // memory determined by profile
    time { (record.n_reads / 300000000).round(2) * 8.hour * params.time_scale_factor }
    module { "${record.tool}/${record.tool_version}" }
    //container "library://singlecell/${record.tool}:${record.tool_version}"

    input:
      val record
    output:
      tuple val(record), path("${record.tool_pubdir}/*"), emit: spaceranger_outputs
      tuple val(record), path("${record.tool_pubdir}/*"), emit: hash_dir
      tuple val(record), path("${record.tool_pubdir}/*{metrics,summary}*", glob: true), emit: metrics

    script:
    main_options = construct_vis_cli_options(record)
    localmem = Math.round(task.memory.toGiga() * 0.95)

    serial_prefix = record.slide.substring(0, 6)
    spatial_dir = "${record.tool_pubdir}/spatial"

    visium_hd_prefixes = ["H1", "SJ", "14072023", "14082023", "26062023", "RD", "UN"]
    is_visium_hd = visium_hd_prefixes.collect { serial_prefix.startsWith(it) }.any()
    is_spaceranger3 = record.tool_version[0].toInteger() >= 3
    """
    spaceranger count $main_options --localcores=$task.cpus --localmem=$localmem

    mkdir -p ${record.tool_pubdir}/_files
    mv ${record.output_id}/_* ${record.tool_pubdir}/_files
    mv ${record.output_id}/*.tgz ${record.tool_pubdir}/
    mv ${record.output_id}/outs/* ${record.tool_pubdir}/

    # new 2024-06-10
    # Visium HD doesn't utilize GPR files (uses VLF)
    # Moreover, adding symlinks in the deliverable for users
    if [[ "${is_visium_hd}" == "true" ]]; then
        cd ${record.tool_pubdir}
        ln -s binned_outputs/square_008um/filtered_feature_bc_matrix .
        ln -s binned_outputs/square_008um/filtered_feature_bc_matrix.h5 .
        ln -s binned_outputs/square_008um/raw_feature_bc_matrix .
        ln -s binned_outputs/square_008um/raw_feature_bc_matrix.h5 .
        cd -
    else
        find ${record.output_id}/SPATIAL_RNA_COUNTER_CS/ -type f -name "metrics_summary_json.json" -exec mv {} ${record.tool_pubdir}/summary.json \\;
        find ${record.output_id}/SPATIAL_RNA_COUNTER_CS/ -type f -name "alignment_metrics.json" -exec mv {} ${record.tool_pubdir}/spatial/alignment_summary.json \\;
        # new 2024-01-02
        # spaceranger containers wont have gprreader on PATH so build its path manually
        sr_root=\$(which spaceranger | xargs dirname)
        gprreader="\${sr_root}/lib/bin/gprreader"
        \$gprreader fetch ${record.slide} ${spatial_dir} --area=${record.area}

        # the above unfortunately doesn't get the raw GPR file, but just a JSON representation
        # below is to pull the raw GPR

        url_base="http://s3-us-west-2.amazonaws.com/10x.spatial-slides/gpr"
        wget -O "${spatial_dir}/${record.slide}.gpr" "\${url_base}/${serial_prefix}/${record.slide}.gpr"
    fi

    if [[ ${is_spaceranger3} ]]; then
        find ${record.output_id}/SPATIAL_RNA_COUNTER*/ -type f -name "metrics_summary_json.json" -exec mv {} ${record.tool_pubdir}/summary.json \\;
    fi
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
    if [ -f "${record.roi_json}" ]; then
        cp ${record.roi_json} .
    fi
    if [ -f "${record.cyta_image}" ]; then
        cp ${record.cyta_image} .
    fi
    if [ -f "${record.dark_image}" ]; then
        cp ${record.dark_image} .
    fi
    if [ -f "${record.color_image}" ]; then
        cp ${record.color_image} .
    fi
    if [ -f "${record.slide_file}" ]; then
        cp ${record.slide_file} .
    fi
    if [ -f "${record.manual_alignment}" ]; then
        cp ${record.manual_alignment} .
    fi
    if [ -f "${record.raw_image}" ]; then
        cp ${record.raw_image} .
    fi
    md5sum * > hashes_img.md5
    """
}
