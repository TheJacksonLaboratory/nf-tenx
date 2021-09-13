#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { compute_fastq_hashes; compute_processed_hashes } from '../modules/hashes.nf'
include { run_citeseq_count } from '../modules/citeseq_count.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'

workflow CITESEQ {
    take: citeseq_records
    main:
    citeseq_records.view{ record -> "CITESEQ / HTO Record: $record.output_id, $record" }
    compute_fastq_hashes(citeseq_records)
    FASTQC(citeseq_records)
    MULTIQC(FASTQC.out.fastqc_results)

    run_citeseq_count(citeseq_records)
    compute_processed_hashes(run_citeseq_count.out.hash_dir)
}
