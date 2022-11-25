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
include { TENX_MULTI } from './tenx_multi'
include { CITESEQ } from './citeseq'
include { INPUT_CHECK } from './input_preprocess'


if (params.samplesheet) { 
    ch_samplesheet = channel.fromPath(params.samplesheet)
} else { exit 1, "Samplesheet not specificed!"}


workflow TENX {
    ch_versions = Channel.empty()

    INPUT_CHECK(
        ch_samplesheet
    )
    ch_versions.mix(INPUT_CHECK.out.versions)
    
    records = INPUT_CHECK.out.main_records
        .branch{
            arc: it.tool == "cellranger-arc"
            atac: it.tool == "cellranger-atac"
            citeseq: it.tool == "citeseq-count"
            gex: (it.tool == "cellranger") && (it.command in ["count", null]) 
            vdj: (it.tool == "cellranger") && (it.command == "vdj")
            multi: (it.tool == "cellranger") && (it.command == "multi")
            visium: it.tool == "spaceranger"
            gex_aggr: (it.tool == "cellranger") && (it.command == "aggr")
            atac_aggr: (it.tool == "cellranger-atac") && (it.command == "aggr")
            arc_aggr: (it.tool == "cellranger-arc") && (it.command == "aggr")
        }

    
    TENX_GEX(records.gex)
    TENX_ATAC(records.atac)
    TENX_ARC(records.arc)
    TENX_VISIUM(records.visium)
    TENX_VDJ(records.vdj)
    TENX_MULTI(records.multi)
    CITESEQ(records.citeseq)
}

//workflow.onComplete {
//    if (params.email || params.email_on_fail) {
//        EmailTemplate.email(workflow, params, summary_params, projectDir, log)
//    }
//    //NfcoreTemplate.summary(workflow, params, log, fail_percent_mapped, pass_percent_mapped)
//}
