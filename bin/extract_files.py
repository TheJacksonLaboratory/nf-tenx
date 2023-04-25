#!/usr/bin/env python3

from pathlib import Path

from arg_utils import parse_cl_paths

args = parse_cl_paths(
    (
        '--ref-path',
        '-r',
        'The reference genome path, as supplied by the "reference_path" field in the samplesheet.',
    ),
    (
        '--annotated_genes_dir',
        '-a',
        'the assets directory containing parquet files, each corresponding to a genome. Each parquet file contains all of the genes of the genome and their gene-types. See the HTML summary document for more information.',
    ),
)
annot_paths = (
    path
    for path in args.gene_annotations_dir.iterdir()
    if path.stem in str(args.ref_genome_dir)
)

with Path('found_genomes.txt').open('w') as f:
    for actual_path in annot_paths:
        f.write(f'{actual_path.stem}\n')
        (output_dir / Path(actual_path.name)).symlink_to(actual_path)

for path in Path.cwd().glob('**/*'):
    print(str(path))
output_dirs = {'cellranger': Path('cellranger')}
output_dir = Path('extracted_files')
output_dir.mkdir()

# Look in each directory for the relevant files and ensure there is only one.
# If so, symlink to it to expose it to the rest of the pipeline
dir_to_pattern = {
    Path.cwd(): '*filtered*matrix*.h5*',
    args.ref_genome_dir: '*genes*.gtf',
}
for direc, pattern in dir_to_pattern.items():
    paths = list(direc.glob(f'**/{pattern}'))
    if len(paths) > 1:
        raise ValueError(
            f'{", ".join(str(path) for path in paths)} were found to match the glob pattern {pattern}. This pipeline '
            f'requires exactly one file of this type to match this pattern.'
        )
    actual_path = paths[0]
    (output_dir / Path(actual_path.name)).symlink_to(actual_path)
