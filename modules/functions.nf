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
  it.collect { it.value != null ? /$it.key="$it.value"/ : /$it.key/ } join " "
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


def construct_cellplex_library_csv_content(record) {
  major_version = record.tool_version[0].toInteger()
  if (major_version >= 9) {
    bam_row = "create-bam,${!record.no_bam}"
  } else {
    bam_row = "no-bam,${record.no_bam ?: ''}"
  }

  rows = [
    "[gene-expression]",
    "reference,${record.reference_path}",
    "expect-cells,${record.n_cells}",
    "include-introns,${record.is_nuclei}",
    bam_row,
    "[libraries]", 
    "fastq_id,fastqs,feature_types"
  ]

  nlibs = record.libraries.size()
  nfqps = record.fastq_paths.size()
  ntypes = record.library_types.size()
  if ([nlibs, nfqps, ntypes].toSet().size() > 1) {
    throw new Exception("Currently can't do unequal #libs, #fastq paths, #types")
  }
  for (i in 0..<nlibs) {
    rows.add("${record.prefixes[i]},${record.fastq_paths[i]},${record.library_types[i]}")
  }

  rows << "[samples]"
  rows << "sample_id,cmo_ids,description"
  record.design.each { cmos, info -> rows.add("${info.name},${cmos},${info.description}") }

  return(rows.join("\n"))
}


def construct_flex_library_csv_content(record) {
  major_version = record.tool_version[0].toInteger()
  if (major_version >= 9) {
    bam_row = "create-bam,${!record.no_bam}"
  } else {
    bam_row = "no-bam,${record.no_bam ?: ''}"
  }

  probeset = file("${params.probe_dir}/${record.probe_set}", checkIfExists: true)
  rows = [
    "[gene-expression]",
    "reference,${record.reference_path}",
    "probe-set,${probeset.toString()}",
    bam_row,
    "[libraries]", 
    "fastq_id,fastqs,feature_types"
  ]

  nlibs = record.libraries.size()
  nfqps = record.fastq_paths.size()
  ntypes = record.library_types.size()
  if ([nlibs, nfqps, ntypes].toSet().size() > 1) {
    throw new Exception("Currently can't do unequal #libs, #fastq paths, #types")
  }
  for (i in 0..<nlibs) {
    rows.add("${record.prefixes[i]},${record.fastq_paths[i]},${record.library_types[i]}")
  }

  rows << "[samples]"
  rows << "sample_id,probe_barcode_ids,description"
  record.design.each { pid, info -> rows.add("${info.name},${pid},${info.description}") }

  return(rows.join("\n"))
}


def create_feature_reference(record) {
    content = ["id,name,read,pattern,sequence,feature_type"]
    feature_type = record.get("feature_type")
    params.tag_list.eachLine { line ->
        def (tag_type, tag_id, read_num, offset5p, tag_sequence, tag_name, pattern) = line.split(",")
        //if (tag_type in record["library_types"] && tag_id in record["tags"]) {
        if (tag_id in record["tags"]) {
            content.add("${tag_id},${tag_name},${read_num},${pattern},${tag_sequence},${feature_type}")
        }
    }
    return content.join("\n")
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
