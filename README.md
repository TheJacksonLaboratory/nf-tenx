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
nextflow run -profile elion TheJacksonLaboratory/nf-tenx --samplesheet path/to/samplesheet.yml
```
### Optional Parameters for single cell gene expression
The following parameters can optionally be passed in for single cell gene expression data:
- `--calc-rna-velo`: using this parameter will run `velocyto run10x --logic=ObservedSpanning10X` on the data. Other `velocyto logic` options will be supported in the future. The resulting `.loom` file has 3 layers, corresponding to `spliced`, `unspliced`, and `ambiguously_spliced` RNA molecules. The pipeline automatically adds them as layers to the resulting `AnnData` (`.h5ad`) object and as assays to the resulting `Seurat` (`.rds`) object. To run with this option:
  ```{bash}
  nextflow run -profile elion TheJacksonLaboratory/nf-tenx --samplesheet path/to/samplesheet.yml --calc-rna-velo
  ```
- `--annots_dir`: for humans and mice (`GRCh37`, `GRCh38`, `GRCm38`, and `GRCm39`) the pipeline automatically [annotates](bin/annotate.py) mitochondrial, ribosomal, hemoglobin, sex-linked, cell-cycle, and stress response genes using `pandas.DataFrame`s in [the `ref_annotations` directory](assets/ref_annotations/), which were created according to this [methodology](#gene-annotation-methodology). This also works with "barnyard" experiments, where your sample is a mix of human and mouse. Other species will be supported in the future. If you'd like to have your data annotated using your own annotations (whether it is a different species or human/mouse), you can pass in your own annotation directory. To run with your own annotations:
  ```{bash}
  nextflow run -profile elion TheJacksonLaboratory/nf-tenx --samplesheet path/to/samplesheet.yml --annots_dir path/to/ref_annotations_directory/
  ``` 
  This directory must contain one or more files that can be read into `pandas.DataFrame` objects. The files must meet the following requirements in order to work:
  - The files must be one of the following file-types: `.parquet`, `.csv`, `.feather`, or `.pickle`. More file types will be supported in the future.
  - When read by `pandas` into a `pandas.DataFrame`, there must be a column of Ensembl gene IDs called "`ensembl_gene_id`". There must also be one or more `boolean` columns for each type of gene you'd like annotated. For example, in the default `.parquet` files that come with this pipeline, there is a column called `mitochondrial` with `boolean` values for all the genes in the genome.
  - The file name must contain a genome (such as `GRCh38`) and an accepted file extension. It may contain other strings (for example, one of the pipeline's defaults is called [`GRCh38.p13.parquet`](assets/ref_annotations/GRCh38.p13.parquet)). Note that only file names containing genomes corresponding to the sample's reference genomes will be used as gene annotations. For example, if your sample's reference genomes are `GRCh38` and `mRatBN7`, and the directory you pass in contains:
    - `GRCh38.feather`
    - `GRCh37.csv`
    - `GRCh38_other_genes.p13.pickle`
    - `GRCm39.parquet`
    - `extraneous_dataframe.pickle`
    - `mRatBN7.2.csv`

    the files `GRCh38.feather`, `GRCh38_other_genes.p13.pickle`, and `mRatBN7.2.csv` will be used to annotate your data based on the correspondence between the gene IDs in your data and the `ensembl_gene_id` column in the files. Again, this will work with any number of files and species.

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
## Gene Annotation Methodology
- Mitochondrial Genes: The [ENSEMBL database](https://ensembl.org) provides the chromosome name for all human and mouse genes. Genes on the mitochondrial chromosome (as identified by the ENSEMBL database) are labeled as `mitochondrial`.
- Hemoglobin Genes: The ENSEMBL database provides descriptions for many genes. Genes containing the world 'hemoglobin' in the description were labeled as `hemoglobin`.
- Sex-linked Genes: Genes on the Y-chromosome in the ENSEMBL database were labeled as `sex_linked`, along with the X-inactivation gene, known as *Xist* in both mice and humans.
- Ribosomal Genes: Genes with a name that starts with the following strings were annotated as `ribosomal`: `Mrpl`, `Mrps`, `Rpl`, `Rps`.
- Cell Cycle Genes: For humans, genes from [this spreadsheet](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6788831/bin/mjy063_supplementary_table_s3.xlsx), a part of [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6788831/) were labeled as `cell_cycle`. Note that there is a column of 'manual annotations' - genes that have a value of 'uncharacterised' in this column were ommitted from the annotation. For mice, all genes in the first two columns of [this spreadsheet](https://static-content.springer.com/esm/art%3A10.1038%2Fnature20123/MediaObjects/41586_2016_BFnature20123_MOESM100_ESM.xlsx) from [this paper](https://www.nature.com/articles/nature20123#Sec23) were used.
- Stress Response Genes: For humans, the genes in [this spreadseet](https://github.com/kieranrcampbell/scrnaseq-digestion-paper/blob/master/data/deliverables/coregene_df-FALSE-v3.csv) from [this paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1830-0) were labeled as `stress_response`. For mice, the genes resulting from [this search](https://www.informatics.jax.org/go/term/GO:0033554) were used.
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
