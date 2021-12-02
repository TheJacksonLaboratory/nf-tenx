#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { COMPUTE_FASTQ_HASHES; COMPUTE_PROCESSED_HASHES } from '../modules/hashes.nf'
include { CELLRANGER_VDJ } from '../modules/cellranger_vdj.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { DUMP_METADATA } from '../modules/metadata.nf'


workflow TENX_VDJ {
    take: vdj_records
    main:
    COMPUTE_FASTQ_HASHES(vdj_records)
    FASTQC(vdj_records)
    MULTIQC(FASTQC.out.fastqc_results)

    CELLRANGER_VDJ(vdj_records)
    COMPUTE_PROCESSED_HASHES(CELLRANGER_VDJ.out.hash_dir)

    metadata_input = COMPUTE_FASTQ_HASHES.out.input_hashes
        .join(COMPUTE_PROCESSED_HASHES.out.output_hashes, remainder:true)
        .join(CELLRANGER_VDJ.out.metrics, remainder:true)
    DUMP_METADATA(metadata_input)
}
