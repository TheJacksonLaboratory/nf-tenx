#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { compute_fastq_hashes; compute_processed_hashes } from '../modules/hashes.nf'
include { run_cellranger_count } from '../modules/cellranger_gex.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'

workflow TENX_GEX {
    take: gex_records
    main:
    gex_records.view{ record -> "Record: $record.output_id, $record" }
    compute_fastq_hashes(gex_records)
    FASTQC(gex_records)
    MULTIQC(FASTQC.out.fastqc_results)

    run_cellranger_count(gex_records)
    compute_processed_hashes(run_cellranger_count.out.hash_dir)
}
