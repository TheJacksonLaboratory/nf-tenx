#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process GET_SPECIES_AND_GENES_GTF {
    tag "$record.output_id"
    executor 'local'

    input:
    val record

    output:
    tuple val(record), path('*')
    path '*species*.txt', emit: species_path
    path '*genes*.gtf', emit: genes_gtf

    script:
    """
    get_species_and_genes_gtf.py --ref_genome_dir=${record.reference_path}
    """
    }