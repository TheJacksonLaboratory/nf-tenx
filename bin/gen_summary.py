#!/usr/bin/env python3

from csv import DictReader
from pathlib import Path

from arg_utils import parse_cl_paths
from directory_tree import display_tree
from jinja2 import Environment, FileSystemLoader


def gen_summary(pubdir: Path, summary_dir: Path, plots_dir: Path) -> None:
    # Get the plots and summary files
    plot_paths = tuple(plots_dir.iterdir())
    summary_files = tuple(summary_dir.iterdir())

    # The CSVs contain information, so get those separately
    extracted_info = {}
    csvs = (path for path in summary_files if path.suffix == '.csv')
    for csv in csvs:
        with csv.open() as f:
            extracted_info[csv.name] = tuple(DictReader(f))

    # Get the flowchart path from the summary files
    flowchart_path = next(path for path in summary_files if path.suffix == '.svg')

    # Use the name of the file and symlink to the actual flowchart to expose to HTML summary
    Path(flowchart_path.name).symlink_to(flowchart_path)

    # Create a directory tree to pass into HTML template
    dir_tree = display_tree(str(pubdir), string_rep=True)

    # Define a list of page sections to pass into the template
    page_sections = [
        'Overview',
        'Plots',
        'Outputs',
        'Annotation Strategy',
        'Tools Used',
    ]

    templates_dir = next(path for path in summary_files if 'templates' in path.name)
    # Create template and HTML and write to a file
    env = Environment(loader=FileSystemLoader(templates_dir))
    template = env.get_template((next(path.name for path in templates_dir.iterdir())))
    content = template.render(
        pipeline_info=extracted_info['overview.csv'],
        page_sections=page_sections,
        plot_paths=plot_paths,
        flowchart=flowchart_path.name,
        dir_tree=dir_tree,
        annotation_descriptions=extracted_info['annotation_descriptions.csv'],
        tools_used=extracted_info['tools.csv'],
    )
    with Path('summary.html').open(mode='w') as f:
        f.write(content)


if __name__ == '__main__':
    args = parse_cl_paths(
        {'--pubdir': '-p', '--summary_dir': '-s', '--plots_dir': '-l'}
    )
    gen_summary(
        pubdir=args.pubdir, summary_dir=args.summary_dir, plots_dir=args.plots_dir
    )
