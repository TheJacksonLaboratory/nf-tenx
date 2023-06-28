#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process FILTER_AMBIENT_RNA {
    tag "$record.output_id-soupx"
    label 'alternative_count'
    publishDir "${params.pubdir}/${record.output_id}/annotations", pattern: "soupx/*filtered*", mode: "copy"
    container 'library://singlecell/r4.2-soupx-seurat:latest'

    input:
    tuple val(record), path('*', stageAs: "cellranger/*")

    output:
    tuple val(record), path('soupx/*', includeInputs: true), val("soupx")

    script:
    """
    mkdir soupx
    filter_ambient_rna.r --args ${record.tool_version} cellranger soupx
    rm -r cellranger/*filtered*
    mv cellranger/* soupx
    """
}
