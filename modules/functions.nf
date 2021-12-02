#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/
import groovy.json.*
import org.yaml.snakeyaml.Yaml


//
// From nf-core/rnaseq
// Extract name of software tool from process name using $task.process
//
def getSoftwareName(task_process) {
    return task_process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()
}

//
// From nf-core/rnaseq
// Extract name of module from process name using $task.process
//
def getProcessName(task_process) {
    return task_process.tokenize(':')[-1]
}



def join_map_items(it) {
  it.collect { it.value ? /$it.key="$it.value"/ : /$it.key/ } join " "
}


def construct_library_csv_content(record) {
  rows = ["fastqs,sample,library_type"]
  nlibs = record.libraries.size()
  nfqps = record.fastq_paths.size()
  ntypes = record.library_types.size()
  if ([nlibs, nfqps, ntypes].toSet().size() > 1) {
    throw new Exception("Currently can't do unequal #libs, #fastq paths, #types")
  }
  for (i in 0..<nlibs) {
    rows << "${record.fastq_paths[i]},${record.prefixes[i]},${record.library_types[i]}"
  }
  return(rows.join("\n"))
}


def compute_runtime(record) {
    n_reads = record.n_reads
    frac = n_reads.divide(3e8)
    runtimes_per_300m = [
        "cellranger-count": 10,
        "cellranger-vdj": 3,
        "cellranger-atac-count": 18,
        "cellranger-arc-count": 20,
        "spaceranger-count": 10
    ]
    key = "${record.tool}-${record.command}"
    if (!(key in runtimes_per_300m)) {
        return("24h")
    }
    return((runtimes_per_300m[key] * frac).round(1))
}
