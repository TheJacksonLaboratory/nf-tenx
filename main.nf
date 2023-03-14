#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    EMERGENCY PARAM SETUP
========================================================================================
*/

params.pubdir = params.getOrDefault("pubdir", "pubdir")
params.publish_fastqs = params.getOrDefault("publish_fastqs", true)
params.probe_dir = params.getOrDefault(
    "probe_dir", 
    "${workflow.projectDir}/assets/probe_sets"
)

/*
========================================================================================
    INITIALIZE
========================================================================================
*/

WorkflowMain.initialize(workflow, params, log)

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

include { TENX } from './workflows/tenx'


workflow NF_TENX {
    TENX()
}

/*
========================================================================================
    RUN WORKFLOWS
========================================================================================
*/

workflow {
    NF_TENX()
}
