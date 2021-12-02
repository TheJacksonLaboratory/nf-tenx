#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { COMPUTE_FASTQ_HASHES; COMPUTE_PROCESSED_HASHES } from '../modules/hashes.nf'
include { CELLRANGER_ARC_COUNT } from '../modules/cellranger_arc.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { DUMP_METADATA } from '../modules/metadata.nf'


workflow TENX_ARC {
    take: arc_records
    main:
    COMPUTE_FASTQ_HASHES(arc_records)
    FASTQC(arc_records)
    MULTIQC(FASTQC.out.fastqc_results)

    CELLRANGER_ARC_COUNT(arc_records)
    COMPUTE_PROCESSED_HASHES(CELLRANGER_ARC_COUNT.out.hash_dir)

    metadata_input = COMPUTE_FASTQ_HASHES.out.input_hashes
        .join(COMPUTE_PROCESSED_HASHES.out.output_hashes, remainder:true)
        .join(CELLRANGER_ARC_COUNT.out.metrics, remainder:true)
    DUMP_METADATA(metadata_input)
}
