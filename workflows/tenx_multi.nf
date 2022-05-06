#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { COMPUTE_FASTQ_HASHES; COMPUTE_PROCESSED_HASHES } from '../modules/hashes.nf'
include { CELLRANGER_MULTI } from '../modules/cellranger_multi.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { SEQUENCING_SATURATION } from '../modules/saturation.nf'
include { DUMP_METADATA } from '../modules/metadata.nf'


workflow TENX_MULTI {
    take: multi_records
    main:
    COMPUTE_FASTQ_HASHES(multi_records)
    FASTQC(multi_records)
    MULTIQC(FASTQC.out.fastqc_results)

    CELLRANGER_MULTI(multi_records)
    COMPUTE_PROCESSED_HASHES(CELLRANGER_MULTI.out.hash_dir)

    metadata_input = COMPUTE_FASTQ_HASHES.out.input_hashes
        .join(COMPUTE_PROCESSED_HASHES.out.output_hashes, remainder:true)
        .join(CELLRANGER_MULTI.out.metrics, remainder:true)
    DUMP_METADATA(metadata_input)
}
