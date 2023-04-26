#!/usr/bin/env python3

from arg_utils import parse_cl_paths
from argparse import Namespace
from pathlib import Path
import scanpy as sc


def extract_files(args: Namespace) -> str:
    """
    Takes reference genome path, a directory of pre-determined gene-annotations, and a CSV path and prepares/stages
    files for next step of pipeline.

    This function does three main tasks:

    - Check the reference genome path in the samplesheet and get the gene annotation files for the genomes present,
    symlinking them to expose them to the rest of the pipeline. It also writes these found genomes to a text file for
    later reporting.

    - Checks that there is just one genes.gtf file and one filtered feature matrix file. If these conditions are met,
    symlink the former to expose it to the rest of the pipeline for velocyto to use.

    - Checks that the genomes present in the filtered feature matrix file match those of the reference genome path.
    If this condition is met, write a CSV mapping each gene symbol to its genome for soupx to use later.
    """

    # Use generator comprehension to get all the gene annotation files of the genomes present in the reference genome
    # path
    annot_paths = (
        path
        for path in args.annots_dir.iterdir()
        if path.stem in str(args.ref_genome_dir)
    )

    # Iterate over each gene annotation file and write the genome found to a file for later reporting. Also symlink
    # it to expose to the rest of the pipeline
    with Path('found_genomes.txt').open('w') as f:
        for actual_path in annot_paths:
            f.write(f'{actual_path.stem}\n')
            Path(actual_path.name).symlink_to(actual_path)

    # Define a pattern for each desired file type
    filt_mtx_pattern = '*filtered*matrix*.h5'
    genes_gtf_pattern = '*genes*.gtf'

    # Define a mapping between the directory to search for a given pattern and the pattern itself
    dir_to_pattern = {
        filt_mtx_pattern: Path.cwd(),
        genes_gtf_pattern: args.ref_genome_dir
    }

    # Use dict comprehension to get a list of paths in each directory that match the pattern. Each list should have
    # just one element
    extracted_paths = {
        pattern: list(direc.rglob(pattern))
        for pattern, direc in dir_to_pattern.items()
    }

    # If any of the above searching found more than one file, raise an error
    if any(len(paths) > 1 for paths in extracted_paths.values()):
        raise ValueError(
            'This pipeline requires exactly one file to match each of the following glob patterns:\n'
            + '\n'.join(extracted_paths.keys())
            + '\nThe following were found for each pattern, respectively:\n'
            + '\n'.join(extracted_paths.values())
        )

    # Read the filtered feature matrix file into memory and get the unique genomes
    adata = sc.read_10x_h5(extracted_paths[filt_mtx_pattern][0])
    filt_mtx_genomes = adata.var['genome'].unique()

    # Use list comprehension to get those genomes that aren't in the reference genome path
    mismatch = [
        genome for genome in filt_mtx_genomes if genome not in args.ref_genome_dir
    ]

    # If this list has any elements, there was a problem with cellranger count, so raise an error
    if mismatch:
        raise ValueError(
            f'The following genomes were present in adata.var["genome"] when {actual_path.name} was read with '
            f'adata = sc.read_10x_h5({actual_path.name}). These are not in the reference genome path '
            f'{args.ref_genome_dir}:\n' + '\n'.join(mismatch)
        )

    # Soupx will require a CSV mapping genes to genomes, so write the CSV now. Note index_label=False for better
    # importing into R, per the function's documentation
    adata.var[['gene_ids', 'genome']].to_csv(args.genes_genome_csv, index_label=False)

    genes_gtf_path = extracted_paths[genes_gtf_pattern][0]
    genes_gtf_path.symlink_to(genes_gtf_path.name)

    return 'The filtered feature matrix and genes.gtf files were found successfully. The genomes present in the ' \
           'former matched those in the reference genome path.'


if __name__ == '__main__':
    cl_paths = parse_cl_paths(
        (
            '--ref_genome_dir',
            '-g',
            'The reference genome path of the sample, as supplied by the "reference_path" field in the samplesheet.'
        ),
        (
            '--annots_dir',
            '-a',
            'The assets directory containing parquet files, each corresponding to a genome. Each parquet file contains '
            'the genes in the genomes associated with their gene-types. See the HTML summary document for more '
            'information.'
        ),
        (
            '--genes_genome_csv',
            '-c',
            'The path of the outputted CSV file that will contain a mapping between each gene to its genome for use '
            'by soupx.'
        )
    )
    message = extract_files(cl_paths)
