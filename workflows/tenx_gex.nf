#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include { COMPUTE_FASTQ_HASHES; COMPUTE_PROCESSED_HASHES } from '../modules/hashes.nf'
include { CELLRANGER_COUNT } from '../modules/cellranger_gex.nf'
include { FASTQC; MULTIQC } from '../modules/qc.nf'
include { SEQUENCING_SATURATION } from '../modules/saturation.nf'
include { DUMP_METADATA } from '../modules/metadata.nf'
include { EXTRACT_FILES } from '../modules/extract.nf'
include { VELOCYTO } from '../modules/velocyto.nf'
include { ANNOTATE_WITH_VELO; ANNOTATE_NO_VELO } from '../modules/annotate.nf'
include { PREPARE_SEURAT; CONVERT_TO_SEURAT} from '../modules/seurat.nf'
include { GEN_PLOTS; GEN_SUMMARY } from '../modules/gen_info.nf'
include { FILTER_AMBIENT_RNA } from '../modules/ambient_rna.nf'

workflow TENX_GEX {
    take: gex_records
    main:
    COMPUTE_FASTQ_HASHES(gex_records)
    FASTQC(gex_records)
    MULTIQC(FASTQC.out.fastqc_results)

    CELLRANGER_COUNT(gex_records)
    EXTRACT_FILES(CELLRANGER_COUNT.out.cellranger_outputs, 'genes_by_genome.csv', 'velocyto_genes.gtf')

    //FILTER_AMBIENT_RNA(EXTRACT_FILES.out.extracted)

    // matrices_ch = CELLRANGER_COUNT.out.cellranger_outputs.mix(FILTER_AMBIENT_RNA.out)

    // if (params.calc_rna_velo) {
    //     summary_dir = params.annotation_info_dir / 'with_rna_velo'
    //     velocyto_input = CELLRANGER_COUNT.out.cellranger_outputs.join(EXTRACT_FILES.out, remainder: true)

    //     VELOCYTO(velocyto_input)

    //     annot_input = EXTRACT_FILES.out.join(VELOCYTO.out)
    //     ANNOTATE_WITH_VELO(annot_input)
    //     adata = ANNOTATE_WITH_VELO.out
    // }

    // else {
    //     summary_dir = params.annotation_info_dir / 'no_rna_velo'

    //     ANNOTATE_NO_VELO(EXTRACT_FILES.out)
    //     adata = ANNOTATE_NO_VELO.out
    // }

    // PREPARE_SEURAT(adata)
    // CONVERT_TO_SEURAT(PREPARE_SEURAT.out)
    // GEN_PLOTS(adata)
    // GEN_SUMMARY(GEN_PLOTS.out, summary_dir)

    // SEQUENCING_SATURATION(CELLRANGER_COUNT.out.cellranger_outputs)

    // hash_input = CELLRANGER_COUNT.out.cellranger_outputs
    //     .join(SEQUENCING_SATURATION.out.hash_data, remainder:true)
    //     .map { it -> [ it[0], it[1..-1].flatten()] }

    // COMPUTE_PROCESSED_HASHES(hash_input)

    // metadata_input = COMPUTE_FASTQ_HASHES.out.input_hashes
    //     .join(COMPUTE_PROCESSED_HASHES.out.output_hashes, remainder:true)
    //     .join(CELLRANGER_COUNT.out.metrics, remainder:true)
    // DUMP_METADATA(metadata_input)
}
