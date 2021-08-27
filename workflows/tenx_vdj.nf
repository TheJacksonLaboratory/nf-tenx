#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { compute_fastq_hashes; compute_processed_hashes } from '../modules/hashes.nf'
include { run_cellranger_vdj_count } from '../modules/cellranger_vdj.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { SEQUENCING_SATURATION } from '../modules/saturation.nf'

workflow TENX_VDJ {
    take: vdj_records
    main:
    vdj_records.view{ record -> "Record: $record.output_id, $record" }
    compute_fastq_hashes(vdj_records)
    FASTQC(vdj_records)
    MULTIQC(FASTQC.out.fastqc_results)

    run_cellranger_vdj_count(vdj_records)
    SEQUENCING_SATURATION(run_cellranger_vdj_count.out.cellranger_vdj_outputs)
    compute_processed_hashes(SEQUENCING_SATURATION.out.hash_data)
}
