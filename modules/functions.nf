#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/
import groovy.json.*
import org.yaml.snakeyaml.Yaml


def join_map_items(it) {
  it.collect { /$it.key="$it.value"/ } join " "
}


def load_samplesheet(unused_stdin) {
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
      vdj: it.tool == "cellranger-vdj"
    }

  return(ch)
}


def collect_fastqs(record) {
  fastqs = []
  prefixes = []
  for(fqp in record.fastq_paths) {
    for(lib in record.libraries) {
      file(fqp).eachFile { 
        fastq -> if(fastq =~ /$lib/) { 
           if (!(fastq in fastqs)) { fastqs.add(fastq) }
           // everything before S\d+_L\d+_[IR]\d+_001.fastq.gz
           prefix = (file(fastq).getName() =~ /(.*)_S\d+_L\d+_[IR]\d+_001.fastq.gz/)[0][1]
           if (!(prefix in prefixes)) { prefixes.add(prefix) }
        }
      }
    }
  }
  record["fastqs"] = fastqs
  record["prefixes"] = prefixes
  
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
