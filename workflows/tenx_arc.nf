#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { compute_fastq_hashes; compute_processed_hashes } from '../modules/hashes.nf'
include { run_cellranger_arc_count } from '../modules/cellranger_arc.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'

workflow TENX_ARC {
    take: arc_records
    main:
    arc_records.view{ record -> "Record: $record.output_id, $record" }
    compute_fastq_hashes(arc_records)
    FASTQC(arc_records)
    MULTIQC(FASTQC.out.fastqc_results)

    run_cellranger_arc_count(arc_records)
    compute_processed_hashes(run_cellranger_arc_count.out.hash_dir)
}
