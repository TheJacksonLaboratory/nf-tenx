#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { compute_fastq_hashes; compute_processed_hashes } from '../modules/hashes.nf'
include { run_cellranger_vdj } from '../modules/cellranger_vdj.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { DUMP_METADATA } from '../modules/metadata.nf'


workflow TENX_VDJ {
    take: vdj_records
    main:
    compute_fastq_hashes(vdj_records)
    FASTQC(vdj_records)
    MULTIQC(FASTQC.out.fastqc_results)

    run_cellranger_vdj(vdj_records)
    compute_processed_hashes(run_cellranger_vdj.out.hash_dir)

    metadata_input = compute_fastq_hashes.out.input_hashes
        .join(compute_processed_hashes.out.output_hashes, remainder:true)
        .join(run_cellranger_vdj.out.metrics, remainder:true)
    DUMP_METADATA(metadata_input)
}
