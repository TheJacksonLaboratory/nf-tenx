#!/usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path

# Define expected args and initialize argument parser
expected_args = {'--ref_genome_dir': '-r', '--gene_annotations_dir': '-g'}
parser = ArgumentParser()

# Iterate over each exepcted arg and add to parser, then parse
for long_name, short_name in expected_args.items():
    parser.add_argument(long_name, short_name, required=True, type=Path)
args = parser.parse_args()

output_dir = Path('extracted_files')
output_dir.mkdir()

# Symlink any genomes in both the reference genome directory and the
# gene annotations directory to expose to rest of pipeline easily
annot_paths = (path for path in args.gene_annotations_dir.iterdir()
               if path.stem in str(args.ref_genome_dir))
with (output_dir / Path('found_genomes.txt')).open('w') as f:
    for actual_path in annot_paths:
        f.write(f'{actual_path.stem}\n')
        Path(actual_path.name).symlink_to(actual_path)

# Look in each directory for the relevant files and ensure there is only one. If so, symlink to it to expose it to the rest of the pipeline
dir_to_pattern = {Path.cwd(): '*filtered*matrix*.h5',
                  args.ref_genome_dir: '*genes*.gtf'}
for direc, pattern in dir_to_pattern.items():
    paths = list(direc.glob(f'**/{pattern}'))
    if len(paths) != 1:
        raise ValueError(
            f'{", ".join(str(path) for path in paths)} were found to match the glob pattern {pattern}. This pipeline requires exactly one file of this type to match this pattern.')
    actual_path = paths[0]
    (output_dir / Path(actual_path.name)).symlink_to(actual_path)
