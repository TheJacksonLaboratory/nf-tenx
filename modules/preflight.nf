#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { getProcessName } from './functions'


process CHECK_INPUT {
    tag "$samplesheet"

    time 5.min
    
    container "library://singlecell/minpyyaml:latest"

    input:
    path samplesheet

    output:
    path "*.yml", emit: yml
    path "version.yml", emit: versions

    script:
    """
    check_samplesheet.py $samplesheet

    sleep 20

    cat <<-END_VERSION > version.yml
    ${getProcessName(task.process)}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


process COUNT_READS {
    cpus 1
    memory 500
    input:
    val(record)

    output:
    tuple val(record), path("n_reads")

    script:
    fastqs = record.fastqs.findAll { it =~ /_R1_/ }.join(" ")
    """
    n_lines=\$(zcat $fastqs | wc -l)
    echo \$((n_lines / 4)) > n_reads
    """
}

