#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { COMPUTE_FASTQ_HASHES; COMPUTE_PROCESSED_HASHES } from '../modules/hashes.nf'
include { SPACERANGER_COUNT; IMAGE_PROCESS } from '../modules/spaceranger.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { DUMP_METADATA } from '../modules/metadata.nf'


workflow TENX_VISIUM {
    take: vis_records
    main:
    COMPUTE_FASTQ_HASHES(vis_records)
    FASTQC(vis_records)
    MULTIQC(FASTQC.out.fastqc_results)

    IMAGE_PROCESS(vis_records)
    SPACERANGER_COUNT(vis_records)
    COMPUTE_PROCESSED_HASHES(SPACERANGER_COUNT.out.hash_dir)

    input_hashes = COMPUTE_FASTQ_HASHES.out.input_hashes
        .mix(IMAGE_PROCESS.out.img_hashes)
        .groupTuple(remainder:true)

    output_hashes = COMPUTE_PROCESSED_HASHES.out.output_hashes

    metadata_input = input_hashes
        .join(output_hashes, remainder:true)
        .join(SPACERANGER_COUNT.out.metrics, remainder:true)
    DUMP_METADATA(metadata_input)
}
