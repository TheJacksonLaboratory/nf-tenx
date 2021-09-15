#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { compute_fastq_hashes; compute_processed_hashes } from '../modules/hashes.nf'
include { run_citeseq_count } from '../modules/citeseq_count.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { DUMP_METADATA } from '../modules/metadata.nf'


workflow CITESEQ {
    take: citeseq_records
    main:
    citeseq_records.view{ record -> "CITESEQ / HTO Record: $record.output_id, $record" }
    compute_fastq_hashes(citeseq_records)
    FASTQC(citeseq_records)
    MULTIQC(FASTQC.out.fastqc_results)

    run_citeseq_count(citeseq_records)
    compute_processed_hashes(run_citeseq_count.out.hash_dir)

    metadata_input = compute_fastq_hashes.out.input_hashes
        .join(compute_processed_hashes.out.output_hashes, remainder:true)
        .join(run_citeseq_count.out.metrics, remainder:true)
    DUMP_METADATA(metadata_input)
}
