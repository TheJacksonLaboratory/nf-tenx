#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { count_reads } from './functions.nf'
 

process COUNT_READS {
    tag "$record.output_id"
    time '15m'

    input:
      val record
    output:
      val(new_record)

    script:
    new_record = count_reads(record)
}
