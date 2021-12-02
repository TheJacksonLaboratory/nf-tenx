#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { COMPUTE_FASTQ_HASHES; COMPUTE_PROCESSED_HASHES } from '../modules/hashes.nf'
include { SPACERANGER_COUNT } from '../modules/spaceranger.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { DUMP_METADATA } from '../modules/metadata.nf'


workflow TENX_VISIUM {
    take: vis_records
    main:
    COMPUTE_FASTQ_HASHES(vis_records)
    FASTQC(vis_records)
    MULTIQC(FASTQC.out.fastqc_results)

    SPACERANGER_COUNT(vis_records)
    COMPUTE_PROCESSED_HASHES(SPACERANGER_COUNT.out.hash_dir)

    metadata_input = COMPUTE_FASTQ_HASHES.out.input_hashes
        .join(COMPUTE_PROCESSED_HASHES.out.output_hashes, remainder:true)
        .join(SPACERANGER_COUNT.out.metrics, remainder:true)
    DUMP_METADATA(metadata_input)
}
