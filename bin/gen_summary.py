#!/usr/bin/env python3

from jinja2 import Environment, FileSystemLoader
from csv import DictReader
from pathlib import Path
from argparse import ArgumentParser
from directory_tree import display_tree

cwd = Path.cwd()
plot_paths = list(cwd.glob('plots/*'))

path_patterns = ('*overview*.csv',
                 '*annotation*descriptions*.csv',
                 '*tools*.txt',
                 '*summary*flowchart*.svg',
                 '*html*templates*')
paths = {pattern.replace('*', ''): next(Path.cwd().glob(pattern)) for pattern in path_patterns}

extracted_info = {}
for key, path in (k, p for k, p in paths.items() if '.' in k and 'svg' not in k):
    with path.open() as f:
        extracted_info[key] = list(DictReader(f))

with paths['summaryflowchart.svg'].open() as f:
    flowchart = f.read()

dir_tree = display_tree('pubdir', string_rep=True)

# Define a list of page sections to pass into the template
page_sections = ['Overview', 'Plots', 'Outputs', 'Annotation Strategy', 'Tools Used']

# Create template and HTML and write to a file
env = Environment(loader=FileSystemLoader(paths['htmltemplates']))
template = env.get_template('summary_template.html')
content = template.render(pipeline_info=extracted_info['overview.csv'], page_sections=page_sections,
                          plot_paths=plot_paths,
                          flowchart=flowchart, dir_tree=dir_tree,
                          annotation_descriptions=extracted_info['annotationdescriptions.csv'],
                          tools_used=extracted_info['tools.txt'])
with open('summary.html', mode='w') as f:
    f.write(content)
