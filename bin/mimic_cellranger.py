#!/usr/bin/env python3

from glob import glob
from pathlib import Path

from arg_utils import parse_cl
from shutil import copy

def mimic_cellranger(soupx_dir: Path) -> None:
    """
    Takes files in the cwd, which should be cellranger's output, and
    copies them all into the soupx directory (excluding the filtered
    matrix because soupx generates its own).
    """
    # Define the glob pattern for filtered feature matrices
    filt_ptrn = '*filtered*'

    # The current directory contains cellranger outputs. Look for
    # anything matching the above-defined pattern so long as it's not in
    # the soupx directory, which is a sub-directory of the current working
    # directory
    cr_filt_mtx_paths = {
        Path(p)
        for p in glob(f'**/{filt_ptrn}', recursive=True)
        if soupx_dir.name not in p
    }

    # Iterate over the filtered matrix paths in the cellranger dir, and for
    # each, get the equivalent soupx matrix path. Then rename it to
    for cr_mtx_path in cr_filt_mtx_paths:
        soupx_equivs = (Path(p) for p in glob(f'{soupx_dir}/**/{filt_ptrn}', recursive=True))
        soupx_equivs = tuple(p for p in soupx_equivs if p.is_file() == cr_mtx_path.is_file())

        if len(soupx_equivs) > 1:
            raise ValueError(
                f'Soupx should have only one file with extension {cr_mtx_path.suffix}, but the following were found:\n'
                + '\n'.join(str(p) for p in soupx_equivs)
            )
        
        soupx_equivs[0].rename(soupx_dir / cr_mtx_path.name)

    # Get all other paths in cellranger directory, and move them into
    # the soupx directory to the cellranger directory. Note that the soupx
    # directory is actually inside the cellranger directory at the moment,
    # so it must be excluded
    other_cr_paths = set(Path().iterdir()) - cr_filt_mtx_paths - {soupx_dir}
    for cr_path in other_cr_paths:
        cr_path.rename(soupx_dir / cr_path.name)


if __name__ == '__main__':
    args = parse_cl(('--soupx_dir', '-s', Path))
    mimic_cellranger(args.soupx_dir)
