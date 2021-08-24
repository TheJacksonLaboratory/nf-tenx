#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { compute_fastq_hashes; compute_processed_hashes } from '../modules/hashes.nf'
include { run_cellranger_atac_count } from '../modules/cellranger_atac.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'

workflow TENX_ATAC {
    take: atac_records
    main:
    atac_records.view{ record -> "Record: $record.output_id, $record" }
    compute_fastq_hashes(atac_records)
    FASTQC(atac_records)
    MULTIQC(FASTQC.outs.fastq_results)

    run_cellranger_atac_count(atac_records)
    compute_processed_hashes(run_cellranger_atac_count.out.hash_dir)
}

