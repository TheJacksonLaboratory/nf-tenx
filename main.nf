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
params.assets_dir = workflow.projectDir / "assets"
params.probe_dir = params.getOrDefault(
    "probe_dir",
    params.assets_dir / "probe_sets"
)
params.tag_list = params.getOrDefault(
    "tag_list",
    file(params.assets_dir / "tags.csv", checkIfExists: true)
)
params.gene_annotations_dir = params.getOrDefault(
    "gene_annotations",
    params.assets_dir / "gene_annotations"
)
params.annotation_info_dir = params.getOrDefault(
    "annotation_info",
    params.assets_dir / "annotation_info"
)
params.calc_rna_velo = ""


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
