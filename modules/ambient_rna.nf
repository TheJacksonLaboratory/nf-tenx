#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process FILTER_AMBIENT_RNA {
    tag "$record.output_id"
    publishDir params.pubdir, pattern: "*soupx*", mode: "copy"
    label 'alternative_count'
    conda 'soupx_env'

    input:
    tuple val(record), path('*')

    output:
    tuple val(record), path('*', includeInputs: true)

    script:
    """
    filter_ambient_rna.r
    rm soupx_filt_mat/genes.tsv
    gzip soupx_filt_mat/*
    cp *filtered*feature*matrix*/features.tsv.gz soupx_filt_mat
    unlink *filtered*feature*matrix*
    """
}

process CONVERT_TO_H5AD {
    tag "$record.output_id"
    publishDir params.pubdir, pattern: "*soupx*", mode: "copy"

    input:
    tuple val(record), path('*')

    output:
    tuple val(record), path('*', includeInputs: true)

    script:
    """
    python3 -c 'import scanpy as sc; adata = sc.read_10x_mtx("soupx_filt_mat"); adata.write("soupx_filtered_feature_matrix.h5ad")'
    """
}