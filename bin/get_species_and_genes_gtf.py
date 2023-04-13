#!/usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path
from re import findall

parser = ArgumentParser()
parser.add_argument('--ref_genome_dir', '-g', required=True, type=Path)
args = parser.parse_args()
ref_genome_dir = args.ref_genome_dir

# Create a mapping between genome name and species
genome_to_species = {
    'mm': 'mus_musculus',
    'GRCm': 'mus_musculus',
    'hg': 'homo_sapiens',
    'GRCh': 'homo_sapiens',
}

# Create the regex for findall()
genome_pattern = r'|'.join(genome_to_species.keys())

# Get the genomes in ref_genome_dir and get a list of species from them
found_genomes = findall(pattern=genome_pattern, string=str(ref_genome_dir))
species = [
    species for genome, species in genome_to_species.items() if genome in found_genomes
]

# If no genome was recognized, raise an error
if not species:
    raise ValueError(
        f'No species could be inferred from the reference genome directory {args.ref_genome_dir}'
    )

# Convert species into a newline-delimited string for printing and writing
species = '\n'.join(species)
print(
    f'The following species were inferred from the the reference genome path {ref_genome_dir}:\n{species}'
)
species_path = Path('species.txt')
with species_path.open(mode='w') as f:
    f.write(species)

# Now get the genes.gtf file for velocyto, verifying that there is only one
file_pattern = '*genes*.gtf'
paths = ref_genome_dir.glob(f'**/{file_pattern}')
paths = list(paths)
symlink_name = file_pattern.replace('*', '')
if len(paths) != 1:
    raise ValueError(
        f'{len(paths)} {symlink_name} files were found in {ref_genome_dir}. Exactly 1 is required.'
    )

# Expose the genes.gtf path to the pipeline with a symlink for velocyto
actual_path = paths[0]
symlink_path = Path(symlink_name)
symlink_path.symlink_to(actual_path)
