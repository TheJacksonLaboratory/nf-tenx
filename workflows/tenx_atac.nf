#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { compute_fastq_hashes; compute_processed_hashes } from '../modules/hashes.nf'
include { run_cellranger_atac_count } from '../modules/cellranger_atac.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { DUMP_METADATA } from '../modules/metadata.nf'


workflow TENX_ATAC {
    take: atac_records
    main:
    compute_fastq_hashes(atac_records)
    FASTQC(atac_records)
    MULTIQC(FASTQC.out.fastqc_results)

    run_cellranger_atac_count(atac_records)
    compute_processed_hashes(run_cellranger_atac_count.out.hash_dir)

    metadata_input = compute_fastq_hashes.out.input_hashes
        .join(compute_processed_hashes.out.output_hashes, remainder:true)
        .join(run_cellranger_atac_count.out.metrics, remainder:true)
    DUMP_METADATA(metadata_input)
}

