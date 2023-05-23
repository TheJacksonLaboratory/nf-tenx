#!/usr/bin/env python3

from glob import glob
from pathlib import Path

import scanpy as sc
from arg_utils import parse_cl


def extract_files(
    ref_genome_dir: Path, genome_csv_name: Path, genes_gtf_name: Path
) -> str:
    """
    Takes reference genome path, a directory of pre-determined
    gene-annotations, and a CSV path and prepares/stages files for next
    step of pipeline.

    This function does two main tasks:

    - Checks that there is just one genes.gtf file and one filtered
    feature matrix file. If these conditions are met, symlink the former
    to expose it to the rest of the pipeline for velocyto to use.

    - Write a CSV mapping each gene symbol to its genome for soupx to
    use later.
    """

    # Define a pattern for each desired file type
    filt_mtx_ptrn = '*filtered*matrix*.h5'
    genes_gtf_ptrn = '*genes*.gtf'

    # Map each pattern to the directory it should be found in
    dir_to_pattern = {
        filt_mtx_ptrn: Path(),
        genes_gtf_ptrn: ref_genome_dir,
    }

    # Use dict comprehension to recursively glob for paths in each
    # directory that match the pattern. Each list should have just one
    # element
    extracted_paths = {
        pattern: glob(f'{direc}/**/{pattern}', recursive=True)
        for pattern, direc in dir_to_pattern.items()
    }

    # If any of the patterns matched more than one file, raise error
    if any(len(paths) > 1 for paths in extracted_paths.values()):
        raise ValueError(
            'This pipeline requires exactly one file to match each of the following '
            'glob patterns:\n'
            + '\n'.join(extracted_paths.keys())
            + '\nThe following were found for each pattern, respectively:\n'
            + '\n'.join(str(extracted_paths.values()))
        )

    # Read the filtered feature matrix file into memory and get the
    # unique genomes
    adata = sc.read_10x_h5(extracted_paths[filt_mtx_ptrn][0])

    # Soupx will require a CSV mapping genes to genomes, so write the
    # CSV now. Note index_label=False for better importing into R, per
    # the function's documentation
    adata.var['genome'].to_csv(genome_csv_name, index_label=False)

    # Finally, symlink the genes.gtf file for velocyto to use later
    genes_gtf_name.symlink_to(extracted_paths[genes_gtf_ptrn][0])

    return (
        'The filtered feature matrix and genes.gtf files were found successfully. The '
        f'mapping between genes and genomes were written to {genome_csv_name}'
    )


if __name__ == '__main__':
    args = parse_cl(
        ('--ref_genome_dir', '-r', Path),
        ('--genome_csv_name', '-c', Path),
        ('--genes_gtf_name', '-g', Path),
    )
    message = extract_files(
        ref_genome_dir=args.ref_genome_dir,
        genome_csv_name=args.genome_csv_name,
        genes_gtf_name=args.genes_gtf_name,
    )
    print(message)
