# nf-tenx

All 10x Genomics data processing pipelines in one containerized workflow

## Purpose

This pipeline is intended to be used at The Jackson Laboratory for processing
10x Genomics data of all kinds.  However, with only a few tweaks to container
locations, this pipeline and its associated tools could be utilized anywhere.

## Supported Assays

Fully supported
- 3'/5' Gene Expression
- 5' Immune Profiling
- Chromatin Accessibility
- Multiome (joint GEX & ATAC)
- CITE-seq with Biolegend TotalSeq-A/B/C ADTs
- Multiplexing with BioLegend TotalSeq-A/C HTOs
- Spatial Transcriptomics (Fresh frozen and FFPE)
- Multiplexing with CellPlex

To do:
- Targetted Gene Expression
- Aggregating all assays

## Usage

1. Configure singularity to pull from JAX's Singularity registry:
```{bash}
srun -q <insert cluster-specific QOS> --time=30 --pty bash
module load singularity
singularity remote add jaxreg jaxreg.jax.org
singularity remote use jaxreg
```
All containers are the `singlecell` collection with the format `<tool-name>:<tool-version>`.

2. Fill out a YAML/JSON-based samplesheet (described below).

3. Run with
```{bash}
nextflow run -profile elion TheJacksonLaboratory/nf-tenx -r 0.3.5 --samplesheet path/to/samplesheet.yml
```

## Samplesheet specifics
The samplesheet format expects a root level YAML sequence or JSON array:
```
# Long form
- sample_name: tinygex
  libraries:
    - LIB123
  library_types:
    - Gene Expression
  ...
# Condensed
- { sample_name: another-sample, libraries: [LIB987], library_types: [Gene Expression], ...}
```

### Required keys
These typically depend on the assay, but most assays require a bare minimum of the following:
- `libraries`: sequence of library names
  each library name must be present somewhere in the FASTQ file names
  e.g. AB12345 would be valid for `XYZ_AB12345_LABNAME_S1_L001_R1_001.fastq.gz`
- `library_types`:
  one of the following corresponding to each library speciified above
  [Gene Expression, Chromatin Accesibility, Spatial Gene Expression, TotalSeq-A, etc.]
- `fastq_paths`: sequeunce of library names.
  a directory path for each of the libraries specified above
- `tool`: a software tool
  [cellranger, spaceranger, cellranger-arc, etc.]
- `tool_version`: a semantic version string for the tool specified above
  e.g. 6.1.1
- `command`: tool command to run
  e.g. `count` to run `cellranger count`
- `sample_name`: any alphanumeric string including hyphens.
- `reference_path`: a file path for the reference transcriptome/genome/etc.

Note that
1. `len(libraries) == len(library_types) == len(fastq_paths)`.
2. even if there is only one library, these three keys must be sequences (e.g. `[lib123]`).

Documentation is in progress for assays that require additional keys, like
TotalSeq (to specify which tags/oligos are in use) or Visium to specify image
and probeset paths.  For now, just look in `modules/assay-type.nf` to see what
cli options are utilized.

Check out the example samplesheets under the `examples` directory.

## Tools available on JAXReg
**Updated 2021-12-02**

| Tool | Version | Notes |
|---|---|---|
| CellRanger | 1.0.0 | |
| | 1.2.0 | |
| | 1.3.0 | |
| | 1.3.1 | |
| | 2.0.2 | |
| | 2.1.0 | |
| | 2.1.1 | |
| | 2.2.0 | |
| | 3.0.0 | |
| | 3.0.2 | |
| | 3.1.0 | |
| | 4.0.0 | |
| | 5.0.0 | |
| | 5.0.1 | |
| | 6.0.0 | |
| | 6.0.1 | |
| | 6.0.2 | |
| | 6.1.1 | |
| | 7.0.0 | |
| CellRanger ATAC | 1.1.0 | do not use, barcode doublet issue |
| | 1.2.0 | |
| | 2.0.0 | |
| CellRanger ARC | 1.0.0 | |
| | 1.0.1 | |
| | 2.0.0 | |
| SpaceRanger | 1.0.0 | |
| | 1.1.0 | |
| | 1.2.0 | |
| | 1.2.2 | |
| | 1.3.0 | |
| | 2.0.0 | |

## Testing

This is a work in progress, but we have compiled and dramatically augmented the
small chr21 test datasets made available from 10x Genomics into a testing suite
for all assays.

To download this data (connected to JAX network):
```{bash}
nextflow pull TheJacksonLaboratory/nf-tenx
nextflow run -profile test nf-tenx/example.nf
```

## Contributing

Please help.
