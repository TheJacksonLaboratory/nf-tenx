#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

import groovy.json.*


def serialize_record(record) {
    return new JsonBuilder(record).toString()
}


def serialize_workflow() {
    def wf = [
        "nextflow_version":   "$nextflow.version",
        "nextflow_build":     "$nextflow.build",
        "nextflow_timestamp": "$nextflow.timestamp",
        "scriptId":       "$workflow.scriptId",
        "scriptName":     "$workflow.scriptName",
        "scriptFile":     "$workflow.scriptFile",
        "repository":     "$workflow.repository",
        "commitId":       "$workflow.commitId",
        "revision":       "$workflow.revision",
        "projectDir":     "$workflow.projectDir",
        "launchDir":      "$workflow.launchDir",
        "workDir":        "$workflow.workDir",
        "homeDir":        "$workflow.homeDir",
        "userName":       "$workflow.userName",
        "configFiles":    "$workflow.configFiles",
        "container":      "$workflow.container",
        "containerEngine":"$workflow.containerEngine",
        "commandLine":    "$workflow.commandLine",
        "profile":        "$workflow.profile",
        "runName":        "$workflow.runName",
        "sessionId":      "$workflow.sessionId",
        "resume":         "$workflow.resume",
        "start":          "$workflow.start",
        "manifest":       "$workflow.manifest"
    ]
    return(new JsonBuilder(wf).toString())
}


process DUMP_METADATA {
    publishDir "${params.pubdir}/${record.output_id}", pattern: "*", mode: "move"

    time 5.min

    container "library://singlecell/minpyyaml:latest"

    input:
    tuple val(record), path("input.md5"), path("output.md5"), path("*")

    output:
    path 'pipeline-metadata.json'

    script:
    workflow_json = serialize_workflow()
    record_json = serialize_record(record)
    """
    echo '${workflow_json}' > workflow.json
    echo '${record_json}' > record.json

    create_metadata.py \
        --record record.json \
        --workflow workflow.json \
        --out pipeline-metadata.json \
        --input-checksums input.md5 \
        --output-checksums output.md5 \
        --metrics *{metrics,summary,report}*
    """
}
