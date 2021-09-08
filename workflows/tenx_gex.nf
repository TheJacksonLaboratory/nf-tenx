#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { compute_fastq_hashes; compute_processed_hashes } from '../modules/hashes.nf'
include { run_cellranger_count } from '../modules/cellranger_gex.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { SEQUENCING_SATURATION } from '../modules/saturation.nf'
include { COUNT_READS } from '../modules/count_reads.nf'

workflow TENX_GEX {
    take: gex_records
    main:
    gex_records.view{ record -> "Record: $record.output_id, $record" }
    gex_records = COUNT_READS(gex_records)
    compute_fastq_hashes(gex_records)
    FASTQC(gex_records)
    MULTIQC(FASTQC.out.fastqc_results)

    run_cellranger_count(gex_records)
    SEQUENCING_SATURATION(run_cellranger_count.out.cellranger_outputs)
    compute_processed_hashes(SEQUENCING_SATURATION.out.hash_data)
}
