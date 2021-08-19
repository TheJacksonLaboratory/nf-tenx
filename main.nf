#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

nextflow.enable.dsl=2

params.pubdir = "pubdir"

include { TENX_GEX } from './workflows/tenx_gex'
include { TENX_ARC } from './workflows/tenx_arc'
include { TENX_ATAC } from './workflows/tenx_atac'
include { TENX_VISIUM } from './workflows/tenx_visium'
include { load_samplesheet } from './modules/functions.nf'

workflow {
    records = load_samplesheet()

    TENX_GEX(records.gex)
    TENX_ATAC(records.atac)
    TENX_ARC(records.arc)
    TENX_VISIUM(records.visium)
}
