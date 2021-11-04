#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { compute_fastq_hashes; compute_processed_hashes } from '../modules/hashes.nf'
include { run_spaceranger_count } from '../modules/spaceranger.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { DUMP_METADATA } from '../modules/metadata.nf'


workflow TENX_VISIUM {
    take: vis_records
    main:
    compute_fastq_hashes(vis_records)
    FASTQC(vis_records)
    MULTIQC(FASTQC.out.fastqc_results)

    run_spaceranger_count(vis_records)
    compute_processed_hashes(run_spaceranger_count.out.hash_dir)

    metadata_input = compute_fastq_hashes.out.input_hashes
        .join(compute_processed_hashes.out.output_hashes, remainder:true)
        .join(run_spaceranger_count.out.metrics, remainder:true)
    DUMP_METADATA(metadata_input)
}
