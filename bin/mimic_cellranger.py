#!/usr/bin/env python3

from arg_utils import parse_cl_paths
from pathlib import Path

def mimic_cellranger(soupx_dir: Path) -> None:
    """
    Takes files in the cwd, which should be cellranger's output, and
    copies them all into the soupx directory (excluding the filtered
    matrix because soupx generates its own). 
    """
    # Define the glob pattern for filtered feature matrices
    filt_ptrn = '*filtered*'

    # The cellranger directory (cr_dir) is the current working directory 
    # (which also contains the soupx directory). Get the filtered feature 
    # matrix paths (a directory and as a .h5 file) and store in a set
    cr_dir = Path()
    cr_filt_mtx_paths = set(cr_dir.glob(filt_ptrn))
    
    # Also get the equivalent soupx paths
    soupx_filt_mtx_paths = tuple(soupx_dir.glob(filt_ptrn))

    # Iterate over the filtered matrix paths in the cellranger dir, and for
    # each, get the equivalent soupx matrix path.
    for cr_mtx_path in cr_filt_mtx_paths:
        soupx_equiv = next(
            soupx_path
            for soupx_path in soupx_filt_mtx_paths
            if cr_mtx_path.suffix == soupx_path.suffix
        )
        
        a = soupx_equiv.rename(soupx_dir / cr_mtx_path.name)

    # Get all other paths in cellranger directory, and symlink them from 
    # the soupx directory to the cellranger directory. Note that the soupx
    # directory is actually inside the cellranger directory at the moment,
    # so it must be excluded
    other_cr_paths = set(cr_dir.iterdir()) - cr_filt_mtx_paths - {soupx_dir}
    for cr_path in other_cr_paths:
        cr_path.rename(soupx_dir / cr_path.name)

if __name__ == '__main__':
    args = parse_cl_paths({'--soupx_dir': '-s'})
    mimic_cellranger(args.soupx_dir)
    