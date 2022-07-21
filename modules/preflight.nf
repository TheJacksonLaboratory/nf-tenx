#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { getProcessName } from './functions'


process CHECK_SAMPLESHEET {
    tag "$samplesheet"

    time 5.min
    executor 'local'
    container "library://singlecell/minpyyaml:latest"

    input:
    path samplesheet

    output:
    path "samplesheet.checked.yml", emit: yml
    path "version.yml", emit: versions

    script:
    """
    check_samplesheet.py $samplesheet samplesheet.checked.yml

    cat <<-END_VERSION > version.yml
    ${getProcessName(task.process)}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


process COUNT_READS {
    cpus 1
    memory 500
    executor 'local'

    input:
    val(record)

    output:
    tuple val(record), path("n_reads")

    script:
    fastqs = record.fastqs.findAll { it =~ /_R1_/ }.join(" ")
    """
    arr=()
    for fq in $fastqs; do arr+=(\$(estFqReadCount -s 5 \$fq)); done
    echo "\${arr[@]/%/+}0" | bc > n_reads
    """
}

