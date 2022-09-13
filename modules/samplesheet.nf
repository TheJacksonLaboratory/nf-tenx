#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/
import groovy.json.*
import org.yaml.snakeyaml.Yaml


def collect_fastqs(LinkedHashMap record) {
  def fastqs = []
  def prefixes = []

  lanes = record.get("lanes")
  if(lanes) {
    lanes = lanes.collect{ it -> it.toString().padLeft(3, "0") }
  } else {
    lanes = ["\\d+"]
  }

  // logic here:
  // for each fastq_path directory, loop over each FASTQ.
  // If that matches as library ID or library ID + lane combination, then add it.
  // Additionally, if `use_undetermined` is specified, add the corresponding
  // Undetermined reads as well
  for(fqp in record.fastq_paths) {
    for(lib in record.libraries) {
      for(lane in lanes) {
        file(fqp).eachFile {

          fastq -> if(fastq =~ /${lib}.*L${lane}.*fastq.*/) {
            if (fastq.isHidden()) { return }
            if (!(fastq in fastqs)) { fastqs.add(fastq.toString()) }

            def prefix = (file(fastq).getName() =~ /(.*)_S\d+_(L\d+_)?[IR]\d+_001.fastq.gz/)[0][1]
            if (!(prefix in prefixes)) { prefixes.add(prefix) }

          } else if(record.get("use_undetermined", false) && fastq =~ /Undetermined.*L${lane}.*fastq.*/ ) {
            fastqs.add(fastq.toString())

            prefix = "Undetermined"
            if (!(prefix in prefixes)) { prefixes.add(prefix) }
          } 

        }
      }
    }
  }
  record["fastqs"] = fastqs.unique()
  record["prefixes"] = prefixes.unique()
  
  return(record)
}


def parse_yaml(yaml_path) {
    // yaml_path is likely a sun.nio.fs.UnixPath
    // or a similar abstract Path if on another OS
    // it needs to be resolved to a string for File to read it
    def text = new FileInputStream(new File(yaml_path.toString()))
    def yaml = new Yaml().load(text)
    return(yaml)
}


def read_n_reads(LinkedHashMap record, f) {
    def text = new Scanner(new File(f.toString())).nextLong()
    record["n_reads"] = text
    return(record)
}


def construct_output_id(LinkedHashMap record) {
    unique_libs = record.libraries.unique()
    record["output_id"] = [unique_libs.join("-"), record.sample_name].join("_")
    return(record)
}


def parse_gt_ids(LinkedHashMap record) {
    def gt_ids = []
    def pattern = ~/GT\d{2}[-_]\d{5}/
    record.fastqs.each { it ->
        def match = it =~ pattern
        if(match) {
            gt_ids += match[0..-1]
        }
    }
    record["GT_IDs"] = gt_ids.toSet().toList()
    return(record)
}


def construct_tool_pubdir(LinkedHashMap record) {
    def pubdir = record.tool
    if (record.command && record.command != "count") {
        pubdir += "-${record.command}"
    }
    record["tool_pubdir"] = pubdir
    return(record)
}

