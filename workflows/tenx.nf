#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { TENX_GEX } from './tenx_gex'
include { TENX_ARC } from './tenx_arc'
include { TENX_ATAC } from './tenx_atac'
include { TENX_VISIUM } from './tenx_visium'
include { TENX_VDJ } from './tenx_vdj'
include { CITESEQ } from './citeseq'
include { load_samplesheet } from '../modules/functions'
include { INPUT_CHECK } from './input_check'


if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, "Samplesheet not specificed!"}


workflow TENX {
    ch_versions = Channel.empty()

    INPUT_CHECK(
        ch_samplesheet
    )
    ch_versions.mix(INPUT_CHECK.out.versions)
    
    records = load_samplesheet(INPUT_CHECK.out.yml)
    
    TENX_GEX(records.gex)
    TENX_ATAC(records.atac)
    TENX_ARC(records.arc)
    TENX_VISIUM(records.visium)
    TENX_VDJ(records.vdj)
    CITESEQ(records.citeseq)
}
