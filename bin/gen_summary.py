#!/usr/bin/env python3

from csv import DictReader
from glob import glob
from pathlib import Path

import pandas as pd
import scanpy as sc
from arg_utils import parse_cl
from directory_tree import display_tree
from jinja2 import Environment, FileSystemLoader
from numpy import nan


def gen_summary(pubdir: Path, summary_dir: Path) -> None:
    # Get the plots and summary files
    plot_paths = glob('**/*.svg', recursive=True)

    # The CSVs contain information, so get those separately
    extracted_info = {}
    for csv in glob(f'{summary_dir}/**/*.csv', recursive=True):
        p = Path(csv)
        with p.open() as f:
            extracted_info[p.name] = tuple(DictReader(f))

    # Get the flowchart path from the summary files
    flowchart_path = Path(glob(f'{summary_dir}/**/*.svg', recursive=True)[0])

    # Use the name of the file and symlink to the actual flowchart to expose to HTML summary
    Path(flowchart_path.name).symlink_to(flowchart_path)

    # Create a directory tree for the pre-analysis steps
    dir_trees = [display_tree(direc, string_rep=True) for direc in pubdir.iterdir()]  # type: ignore

    # Read AnnData object and generate a dict of the special genes
    adata = sc.read_h5ad(glob('**/*.h5ad', recursive=True)[0])
    gene_types = {
        col.title(): adata.var.loc[adata.var[col], :].index.to_list()
        for col in adata.var.columns
        if adata.var[col].dtype == bool
    }

    classes = 'table table-bordered table-light table-striped table-hover align-middle border border-info'
    # Get the longest list of genes in the dict so as to fill it in order
    # to make a pd.DataFrame out of it
    try:
        length = max(len(l) for l in gene_types.values())

        # Fill the lists with np.nan's so they're all the same length
        for gene_list in gene_types.values():
            gene_list.extend([nan] * (length - len(gene_list)))

        # Convert the dict to a dataframe, define CSS class, and convert to
        # HTML
        gene_type_df = pd.DataFrame.from_dict(gene_types)
        gene_annot_table = gene_type_df.to_html(
            index=False, na_rep='', classes=classes, justify='unset'
        )
    except:
        gene_type_df = pd.DataFrame.from_dict(gene_types)
        gene_annot_table = gene_type_df.to_html(
            index=False, na_rep='', classes=classes, justify='unset'
        )

    # TODO - GENERATE THE LINK FOR THE CELLRANGER-GENERATED WEB SUMMARY AND FEED IT TO TEMPLATE

    # Define the title and a list of page sections to pass into the
    # template
    title = 'Pre-Analysis Pipeline'
    page_sections = [
        'Overview',
        'Plots',
        'Outputs',
        'Annotation Strategy',
        'Tools Used',
    ]

    templates_dir = Path(glob(f'{summary_dir}/**/*templates*', recursive=True)[0])

    # Create template and HTML and write to a file
    env = Environment(loader=FileSystemLoader(templates_dir))
    template = env.get_template((next(path.name for path in templates_dir.iterdir())))
    content = template.render(
        title=title,
        pipeline_info=extracted_info['overview.csv'],
        page_sections=page_sections,
        plot_paths=plot_paths,
        flowchart=flowchart_path.name,
        dir_trees=dir_trees,
        annotation_descriptions=extracted_info['annotation_descriptions.csv'],
        gene_annot_table=gene_annot_table,
        tools_used=extracted_info['tools.csv'],
    )
    with open(f'{title.replace(" ", "_")}_Summary.html', mode='w') as f:
        f.write(content)


if __name__ == '__main__':
    args = parse_cl(('--pubdir', '-p', Path), ('--summary_dir', '-s', Path))
    gen_summary(pubdir=args.pubdir, summary_dir=args.summary_dir)
