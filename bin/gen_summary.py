#!/usr/bin/env python3

from jinja2 import Environment, FileSystemLoader
from csv import DictReader
from pathlib import Path
from argparse import ArgumentParser
from directory_tree import display_tree

parser = ArgumentParser()
expected_args = {'--pubdir': '-p', '--summary_dir': '-s'}
for long_name, short_name in expected_args.items():
    parser.add_argument(long_name, short_name, required=True, type=Path)
args = parser.parse_args()

plot_paths = [f'plots/{path.name}' for path in Path.cwd().glob('plots/*')]

path_patterns = ('*overview*.csv',
                 '*annotation*descriptions*.csv',
                 '*tools*.txt',
                 '*summary*flowchart*.svg',
                 '*html*templates*')
paths = {pattern.replace('*', ''): next(args.summary_dir.glob(pattern)) for pattern in path_patterns}

extracted_info = {}
for key, path in ((k, p) for k, p in paths.items() if '.' in k and 'svg' not in k):
    with path.open() as f:
        extracted_info[key] = list(DictReader(f))

# Symlink a path to the flowchart to put it in the pubdir, so the HTML summary has access to it
extracted_info['summaryflowchart.svg'] = Path(paths['summaryflowchart.svg'].name)
extracted_info['summaryflowchart.svg'].symlink_to(paths['summaryflowchart.svg'])

dir_tree = display_tree(args.pubdir, string_rep=True)

# Define a list of page sections to pass into the template
page_sections = ['Overview', 'Plots', 'Outputs', 'Annotation Strategy', 'Tools Used']

# Create template and HTML and write to a file
env = Environment(loader=FileSystemLoader(paths['htmltemplates']))
template = env.get_template('summary_template.html')
content = template.render(pipeline_info=extracted_info['overview.csv'], page_sections=page_sections,
                          plot_paths=plot_paths,
                          flowchart=extracted_info['summaryflowchart.svg'], dir_tree=dir_tree,
                          annotation_descriptions=extracted_info['annotationdescriptions.csv'],
                          tools_used=extracted_info['tools.txt'])
with open('summary.html', mode='w') as f:
    f.write(content)
