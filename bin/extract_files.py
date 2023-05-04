#!/usr/bin/env python3

from arg_utils import parse_cl_paths
from pathlib import Path
import scanpy as sc


def extract_files(
    ref_genome_dir: Path, annots_dir: Path, genome_csv_name: Path, genes_gtf_name: Path
) -> str:
    """
    Takes reference genome path, a directory of pre-determined
    gene-annotations, and a CSV path and prepares/stages files for next
    step of pipeline.

    This function does three main tasks:

    - Check the reference genome path in the samplesheet and get the
    gene annotation files for the genomes present, symlinking them to
    expose them to the rest of the pipeline. It also writes these found
    genomes to a text file for later reporting.

    - Checks that there is just one genes.gtf file and one filtered
    feature matrix file. If these conditions are met, symlink the former
    to expose it to the rest of the pipeline for velocyto to use.

    - Checks that the genomes present in the filtered feature matrix
    file match those of the reference genome path. If this condition is
    met, write a CSV mapping each gene symbol to its genome for soupx to
    use later.
    """

    # Each file in the annotations directory, located in assets,
    # corresponds to a genome. The stem of these files are genomes, so
    # take only those files whose stems are found in the reference
    # genome path. Store them in a dict as keys, and symlink to the
    # files while we're at it using dict comprehension

    # Define a pattern for each desired file type
    filt_mtx_pattern = '*filtered*matrix*.h5'
    genes_gtf_pattern = '*genes*.gtf'

    # Map each pattern to the directory it should be found in
    dir_to_pattern = {
        filt_mtx_pattern: Path(),
        genes_gtf_pattern: ref_genome_dir,
    }

    # Use dict comprehension to recursively glob for paths in each
    # directory that match the pattern. Each list should have just one
    # element
    extracted_paths = {
        pattern: tuple(direc.rglob(pattern))
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
    adata = sc.read_10x_h5(extracted_paths[filt_mtx_pattern][0])
    filt_mtx_genomes = adata.var['genome'].drop_duplicates()

    # Use generator comprehension to get the annotation files whose stems
    # (each representing a genome) match the genomes in the series created
    # above
    annot_paths = (
        path
        for path in annots_dir.iterdir()
        if filt_mtx_genomes.str.match(pat=rf'.*{path.name.split(".")[0]}.*', case=False).any()
    )

    # Iterate and symlink each to expose to the rest of the pipeline
    for actual_path in annot_paths:
        Path(actual_path.name).symlink_to(actual_path)

    # Soupx will require a CSV mapping genes to genomes, so write the
    # CSV now. Note index_label=False for better importing into R, per
    # the function's documentation
    adata.var['genome'].to_csv(genome_csv_name, index_label=False)

    # Finally, symlink the genes.gtf file for velocyto to use later
    genes_gtf_name.symlink_to(extracted_paths[genes_gtf_pattern][0])

    return (
        'The filtered feature matrix and genes.gtf files were found successfully. The '
        'genomes present in the former matched those in the reference genome path.'
    )


if __name__ == '__main__':
    args = parse_cl_paths(
        {
            '--ref_genome_dir': '-r',
            '--annots_dir': '-a',
            '--genome_csv_name': '-c',
            '--genes_gtf_name': '-g',
        }
    )
    message = extract_files(
        ref_genome_dir=args.ref_genome_dir,
        annots_dir=args.annots_dir,
        genome_csv_name=args.genome_csv_name,
        genes_gtf_name=args.genes_gtf_name,
    )
    print(message)
