#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

process SCRUBLET {
}


process SOUPX {
    container "docker://irinashchukina/soupx"

}


process VELOCYTO {
    container "library://singlecell/velocyto:0.17.17"
}


process ANNOTATE_MATRIX {
}
