#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { COMPUTE_FASTQ_HASHES; COMPUTE_PROCESSED_HASHES } from '../modules/hashes.nf'
include { CITESEQ_COUNT } from '../modules/citeseq_count.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { DUMP_METADATA } from '../modules/metadata.nf'


workflow CITESEQ {
    take: citeseq_records
    main:
    COMPUTE_FASTQ_HASHES(citeseq_records)
    FASTQC(citeseq_records)
    MULTIQC(FASTQC.out.fastqc_results)

    CITESEQ_COUNT(citeseq_records)
    COMPUTE_PROCESSED_HASHES(CITESEQ_COUNT.out.hash_dir)

    metadata_input = COMPUTE_FASTQ_HASHES.out.input_hashes
        .join(COMPUTE_PROCESSED_HASHES.out.output_hashes, remainder:true)
        .join(CITESEQ_COUNT.out.metrics, remainder:true)
    DUMP_METADATA(metadata_input)
}
