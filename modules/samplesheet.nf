#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/
import groovy.json.*
import org.yaml.snakeyaml.Yaml


def collect_fastqs(record) {
  fastqs = []
  prefixes = []
  for(fqp in record.fastq_paths) {
    for(lib in record.libraries) {
      file(fqp).eachFile { 
        fastq -> if(fastq =~ /$lib/) { 
           if (!(fastq in fastqs)) { fastqs.add(fastq.toString()) }
           // everything before S\d+_L\d+_[IR]\d+_001.fastq.gz
           prefix = (file(fastq).getName() =~ /(.*)_S\d+_(L\d+_)?[IR]\d+_001.fastq.gz/)[0][1]
           if (!(prefix in prefixes)) { prefixes.add(prefix) }
        }
      }
    }
  }
  record["fastqs"] = fastqs
  record["prefixes"] = prefixes
  
  return(record)
}


def load_sample_sheet(samplesheet) {
    text = new FileInputStream(new File(params.samplesheet))
    yaml = new Yaml().load(text)
    return(yaml)
}

def read_n_reads(record, f) {
    text = new Scanner(new File(f.toString())).nextLong()
    record["n_reads"] = text
    return record
}


def construct_output_id(record) {
    record["output_id"] = [record.libraries.join("-"), record.sample_name].join("_")
    return(record)
}


def parse_gt_ids(record) {
    gt_ids = []
    def pattern = ~/GT\d{2}[-_]\d{5}/
    record.fastqs.each { it ->
        def match = it =~ pattern
        if(match) {
            gt_ids += match[0..-1]
        }
    }
    record["GT_IDs"] = gt_ids.toSet().toList()
    return (record)
}


def construct_tool_pubdir(record) {
    pubdir = record.tool
    if (record.command && record.command != "count") {
        pubdir += "-${record.command}"
    }
    record["tool_pubdir"] = pubdir
    return(record)
}

