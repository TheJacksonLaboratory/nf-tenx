# nf-tenx

All 10x Genomics data processing pipelines in one containerized workflow

## Purpose

This pipeline is intended to be used at The Jackson Laboratory for processing
10x Genomics data of all kinds.  However, with only a few tweaks to container
locations, this pipeline and its associated tools could be utilized anywhere.

## Usage

1. Configure singularity to pull from JAX's Singularity registry:
```{bash}
srun -q <insert cluster-specific QOS> --time=30 --pty bash
module load singularity
singularity remote add jaxreg jaxreg.jax.org
singularity remote use jaxreg
```

2. Fill out a YAML/JSON-based samplesheet:
   -   No YAML documents
   -   Use a sequence/list as the root object.
   -   Example samplesheets are found under `tests/`
   -   Required keys depend on the assay but are typically:
       -   `libraries`: sequeunce of library names.
           each library name must be present somewhere in the FASTQ file names
           e.g. AB12345 would be valid for `XYZ_AB12345_LABNAME_S1_L001_R1_001.fastq.gz`
       -   `library_types`: sequeunce of library types.
           one of the following corresponding to each library speciified above
           [Gene Expression, Chromatin Accesibility, Spatial Gene Expression, Multiplex, etc.]
       -   `fastq_paths`: sequeunce of library names.
           a directory path for each of the libraries specified above
       -   `tool`: a 10x Genomics software tool
           [cellranger, spaceranger, cellranger-arc, etc.]
       -   `tool_version`: a semantic version string for the tool specified above
           e.g. 6.1.1
       -   `sample_name`: any alphanumeric string including hyphens.
       -   `reference_path`: a file path for the reference transcriptome/genome/etc.

3. Run with
```{bash}
nextflow run -profile elion,singularity nf-tenx --samplesheet path/to/samplesheet.yml
```

## Tools available
**Updated 2021-08-27**

CellRanger
- 1.0.0
- 1.2.0
- 1.3.0
- 1.3.1
- 2.0.2
- 2.1.0
- 2.1.1
- 2.2.0
- 3.0.0
- 3.0.2
- 3.1.0
- 4.0.0
- 5.0.0
- 5.0.1
- 6.0.0
- 6.0.1
- 6.0.2
- 6.1.1

CellRanger ATAC
- (<1.2.0 versions are not recommended due to a barcode doublet issue)
- 1.2.0
- 2.0.0
  
CellRanger ARC
- 1.0.0
- 1.0.1
- 2.0.0

SpaceRanger
- 1.0.0
- 1.1.0
- 1.2.0
- 1.2.2
- 1.3.0

## Testing

This is a work in progress, but we have compiled and dramatically augmented the
small chr21 test datasets made available from 10x Genomics into a testing suite
for all assays.

To download this data (connected to JAX network):
```{bash}
nextflow pull TheJacksonLaboratory/nf-tenx
nextflow run -profile singularity,test nf-tenx/example.nf
```

## Contributing
