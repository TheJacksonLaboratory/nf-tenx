#!/usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path
from re import findall

# Define expected args and initialize argument parser
expected_args = [('--ref_genome_dir', '-r'), ('--gene_annotations_dir', '-g')]
parser = ArgumentParser()

# Iterate over each exepcted arg and add to parser, then parse
for long_name, short_name in expected_args:
    parser.add_argument(long_name, short_name, required=True, type=Path)
args = parser.parse_args()

# Get all the files in the gene annotations directory, and generate 
# regex from it (since these files have genome names in them)
supported_genomes = args.gene_annotations_dir.glob('*')
genome_pattern = r'|'.join(genome.stem for genome in supported_genomes)

# Use the regex to get genomes out of the reference genome directory
found_genomes = findall(pattern=genome_pattern, string=str(args.ref_genome_dir))

# Symlink the gene-annotated dataframes for the genomes present in 
# reference data directory to expose them to the pipeline
for actual_path in args.gene_annotations_dir.glob(
    f'*[{",".join(found_genomes)}].pickle'
):
    Path(actual_path.name).symlink_to(actual_path)

# Write the found genomes into a text file for web summary
with Path('found_genomes.txt').open('w') as f:
    f.write('\n'.join(found_genomes))

# Now get the genes.gtf file for velocyto. If there is more than one,
# an error will be thrown
for actual_path in args.ref_genome_dir.glob(f'**/*genes*.gtf'):
    Path(actual_path.name).symlink_to(actual_path)
