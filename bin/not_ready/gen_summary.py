#!/usr/bin/env python3

from jinja2 import Environment, FileSystemLoader
from csv import DictReader
from pathlib import Path
from argparse import ArgumentParser
from directory_tree import display_tree

# Create a list of expected command-line arguments and their settings
# for the argument parser
expected_args = [
    {
        "long_name": "--plots_dir",
        "short_name": "-p",
        "required": True,
        "type": Path,
    },
    {
        "long_name": "--summary_dir",
        "short_name": "-s",
        "required": True,
        "type": Path,
    }
]

# Initialize parser and iterate over each expected argument, adding to
# parser, and then parse the arguments
parser = ArgumentParser()
for arg in expected_args:
    parser.add_argument(
        arg["long_name"],
        arg["short_name"],
        type=arg["type"],
        required=arg["required"]
        )
args = parser.parse_args()

# Define a list of page sections to pass into the template
page_sections = ['Overview', 'Plots', 'Outputs', 'Annotation Strategy', 'Tools Used']

# Get the paths for the plots to be passed into the template
plot_paths = list(args.plots_dir.glob('*'))

# Each step in the pipeline has a description stored in a CSV.
# Get these and load them into memory to pass into the template
pipeline_info = []
overview_path = list(args.summary_dir.glob("**/*overview*.csv"))[0]
with overview_path.open() as f:
    reader = DictReader(f)
    for row in reader:
        pipeline_info.append(row)

# This is pretty much the same thing as the above. Can probably combine
# them both into some kind of iterable
annotation_descriptions = []
annotation_description_path = list(args.summary_dir.glob("**/*gene*annotation*descriptions*.csv"))[0]
with annotation_description_path.open() as f:
    reader = DictReader(f)
    for row in reader:
        annotation_descriptions.append(row)

tools_used = []
tools_path = list(args.summary_dir.glob("**/*tools*.txt"))[0]
with tools_path.open() as f:
    for line in f.readlines():
        tools_used.append(line.strip())

# Get the flowchart by globbing for it and then geting the first
# (and only file)
actual_flowchart_path = list(args.summary_dir.glob("**/*pipeline*summary*.svg"))[0]

# Generate a directory tree to show on summary
im_dir_tree = display_tree('../intermediate_outputs', string_rep=True)
final_dir_tree = display_tree('../final_outputs', string_rep=True)

# Symlink to the actual path to expose the file to the HTML template
flowchart_path = Path("pipeline_summary.svg")
try:
    flowchart_path.symlink_to(actual_flowchart_path)
except FileExistsError:
    pass
# Finally, define the directory where the template lives
templates_dir = list(args.summary_dir.glob("**/*templates*"))[0]

# Create template and HTML and write to a file
env = Environment(loader=FileSystemLoader(templates_dir))
template = env.get_template('summary_template.html')
content = template.render(pipeline_info=pipeline_info, page_sections=page_sections, plot_paths=plot_paths, flowchart_path=flowchart_path, im_dir_tree=im_dir_tree, final_dir_tree=final_dir_tree, annotation_descriptions=annotation_descriptions, tools_used=tools_used)
with open('summary.html', mode='w') as f:
    f.write(content)