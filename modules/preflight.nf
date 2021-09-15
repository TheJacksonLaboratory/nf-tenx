#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process CHECK_INPUT {

    time 5.min
    
    container "library://singlecell/minpyyaml:latest"

    input:
    val samplesheet
    output:
    stdout emit: all_good

    script:
    samplesheet = file(samplesheet)
    println samplesheet
    """
    check_samplesheet.py $samplesheet
    """
}
