#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { CHECK_INPUT } from '../modules/preflight'


workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    CHECK_INPUT(samplesheet)

    emit:
    yml = CHECK_INPUT.out.yml
    versions = CHECK_INPUT.out.versions

}
