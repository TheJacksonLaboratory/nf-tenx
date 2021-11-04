#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { compute_fastq_hashes; compute_processed_hashes } from '../modules/hashes.nf'
include { run_cellranger_count } from '../modules/cellranger_gex.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { SEQUENCING_SATURATION } from '../modules/saturation.nf'
include { DUMP_METADATA } from '../modules/metadata.nf'


workflow TENX_GEX {
    take: gex_records
    main:
    compute_fastq_hashes(gex_records)
    FASTQC(gex_records)
    MULTIQC(FASTQC.out.fastqc_results)

    run_cellranger_count(gex_records)
    SEQUENCING_SATURATION(run_cellranger_count.out.cellranger_outputs)
    compute_processed_hashes(SEQUENCING_SATURATION.out.hash_data)

    metadata_input = compute_fastq_hashes.out.input_hashes
        .join(compute_processed_hashes.out.output_hashes, remainder:true)
        .join(run_cellranger_count.out.metrics, remainder:true)
    DUMP_METADATA(metadata_input)
}
