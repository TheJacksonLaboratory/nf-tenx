#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { compute_fastq_hashes; compute_processed_hashes } from '../modules/hashes.nf'
include { run_spaceranger_count } from '../modules/spaceranger.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'

workflow TENX_VISIUM {
    take: vis_records
    main:
    vis_records.view{ record -> "Record: $record.output_id, $record" }
    compute_fastq_hashes(vis_records)
    FASTQC(vis_records)
    MULTIQC(FASTQC.out.fastqc_results)

    run_spaceranger_count(vis_records)
    compute_processed_hashes(run_spaceranger_count.out.hash_dir)
}
