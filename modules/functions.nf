#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/
import groovy.json.*
import org.yaml.snakeyaml.Yaml


def join_map_items(it) {
  it.collect { it.value ? /$it.key="$it.value"/ : /$it.key/ } join " "
}


def load_samplesheet(unused_stdin) {
  text = new FileInputStream(new File(params.samplesheet))
  yaml = new Yaml().load(text)

  ch = Channel.fromList(yaml)
    .map{
      it -> 
        count_reads_approx(
          collect_fastqs(
            construct_tool_pubdir(construct_output_id(it))
          )
        )
      }
    .branch{
      arc: it.tool == "cellranger-arc"
      atac: it.tool == "cellranger-atac"
      hto: it.tool == "citeseq-count"
      gex: (it.tool == "cellranger")  && (it.command == "count") 
      vdj: (it.tool == "cellranger") && (it.command == "vdj")
      multi: (it.tool == "cellranger") && (it.command == "multi")
      visium: it.tool == "spaceranger"
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


def count_reads_approx(record) {
  // I1 reads are between 8-10bp long (depending on single/dual index and gex vs atac)
  // On average, I've computed each 10bp fastq record to be ~9B compressed (at
  // whatever GT's standard compression level is, hard/impossible to tell)
  // Anyway, at 9B/read, we should have a conservative estimate
  // Idea from: https://gitter.im/nextflow-io/nextflow?at=6139658f5739ab2df8e209d1
  lines = 0
  record.fastqs.each {
    it -> if(it =~ /_I1_/) {
      lines += file(it).size()
    }
  }
  record["n_reads"] = lines.intdiv(9)
  return(record)
}


def construct_output_id(record) {
  record["output_id"] = [record.libraries.join("-"), record.sample_name].join("_")
  return(record)
}


def construct_tool_pubdir(record) {
  pubdir = record.tool
  if (record.command && record.command != "count") {
    pubdir += "-${record.command}"
  }
  record["tool_pubdir"] = pubdir
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
