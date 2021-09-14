#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { compute_fastq_hashes; compute_processed_hashes } from '../modules/hashes.nf'
include { run_cellranger_vdj } from '../modules/cellranger_vdj.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'

workflow TENX_VDJ {
    take: vdj_records
    main:
    vdj_records.view{ record -> "VDJ Record: $record.output_id, $record" }
    compute_fastq_hashes(vdj_records)
    FASTQC(vdj_records)
    MULTIQC(FASTQC.out.fastqc_results)

    run_cellranger_vdj(vdj_records)
    compute_processed_hashes(run_cellranger_vdj.out.hash_dir)
}
