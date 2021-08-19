#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/
import groovy.json.*
import org.yaml.snakeyaml.Yaml


def load_samplesheet() {
  text = new FileInputStream(new File(params.samplesheet))
  yaml = new Yaml().load(text)

  ch = Channel.fromList(yaml)
    .map{
      it -> count_reads(
        collect_fastqs(
          construct_output_id(it)
        )
      )}
    .branch{
      gex: it.tool == "cellranger"
      atac: it.tool == "cellranger-atac"
      arc: it.tool == "cellranger-arc"
      visium: it.tool == "spaceranger"
      hto: it.tool == "citeseq"
    }

  return(ch)
}


def collect_fastqs(record) {
  fastqs = []
  prefixes = []
  for(lib in record.libraries) {
    file(record.fastq_path).eachFile { 
      item -> if(item =~ /$lib/) { 
         fastqs.add(item) 
         // everything before S\d+_L\d+_[IR]\d+_001.fastq.gz
         prefixes.add((file(item).getName() =~ /(.*)_S\d+_L\d+_[IR]\d+_001.fastq.gz/)[0][1])
      }
    }
  }
  record["fastqs"] = fastqs
  record["prefixes"] = prefixes.toSet()
  
  return(record)
}

def count_reads(record) {
  lines = 0
  record.fastqs.each {
    it -> if(it =~ /_I1_/) {
      lines += file(it).countFastq()
    }
  }
  record["n_reads"] = lines
  return(record)
}


def construct_output_id(record) {
  record["output_id"] = [record.libraries.join("-"), record.sample_name].join("_")
  return(record)
}
