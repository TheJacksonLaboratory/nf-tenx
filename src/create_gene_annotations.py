#!/usr/bin/env python3

from io import StringIO
from pathlib import Path
from re import search

import pandas as pd
from bioservices import BioMart
from numpy import loadtxt

mart = 'ENSEMBL_MART_ENSEMBL'

# Archived ensembl datasets exist on other hosts, so map the hosts to the
# appropriate datasets. If None, 
ensembl_hosts_dsets = {
    'grch37.ensembl.org': 'hsapiens_gene_ensembl',
    None: 'hsapiens_gene_ensembl',
    'nov2020.archive.ensembl.org': 'mmusculus_gene_ensembl',
    None: 'mmusculus_gene_ensembl'
}

# Define a tuple of attributes of interest
attributes_of_interest = (
    'ensembl_gene_id',
    'external_gene_name',
    'chromosome_name',
    'description',
)

for host, ds in ensembl_hosts_dsets.items():
    # Establish a connection to ENSEMBL host and set return format to CSV
    s = BioMart(host=host)
    s.custom_query(formatter='CSV')

    # Add the dataset to the XML query and the attributes using list
    # comprehension
    s.add_dataset_to_xml(dataset=ds)
    _ = [s.add_attribute_to_xml(att) for att in attributes_of_interest]

    # Generate the XML query, send request, and store the result, convert
    # to file-like buffer for pandas to read
    xml_query = s.get_xml()
    result_string = s.query(xml_query)
    result_buffer = StringIO(result_string)

    # Use buffer to make a dataframe, where genes will be annotated
    genes_df = pd.read_csv(result_buffer, names=attributes_of_interest)

    # Delete unnecesary items to save memory
    del result_buffer, result_string, _

    # Mitochondrial genes are on the mitochondrial chromosome
    genes_df['mitochondrial'] = genes_df['chromosome_name'].str.match('mt', case=False)

    # Assume that if hemoglobin is in the description, then it's a
    # hemoglobin gene. Also, filter out NaNs
    genes_df['hemoglobin'] = (
        genes_df['description'].str.contains('hemoglobin', case=False)
    ) & (genes_df['description'] == genes_df['description'])

    # Sex linked genes are either on the Y-chromosome or are the
    # x-inactivation gene Xist
    genes_df['sex_linked'] = (
        genes_df['chromosome_name'].str.match('y', case=False)
    ) | (genes_df['external_gene_name'].str.contains('xist', case=False))

    # Ribosomal gene names have a certain pattern, so match that pattern
    # Also filter out NaNs because not all genes have external gene names
    genes_df['ribosomal'] = (
        genes_df['external_gene_name'].str.match('m?rp[ls]', case=False)
    ) & (genes_df['external_gene_name'] == genes_df['external_gene_name'])

    # Cell cycle and stress response genes are taken from pre-built files.
    # Get the paths for each inside the dataset's directory using dict
    # comprehension and globbing
    cc_sr_paths = {
        gene_type.replace('*', '_'): next(
            Path().glob(f'*{ds.removesuffix("_gene_ensembl")}*/*{gene_type}*')
        )
        for gene_type in ('cell*cycle', 'stress*response')
    }

    # Define the reader for each of the three file types present
    readers = {'.csv': pd.read_csv, '.xlsx': pd.read_excel, '.txt': loadtxt}

    # Use dict comprehension to assign a dataframe or numpy array to each
    # gene type in a dict
    cc_sr_genes = {
        gene_type: readers[path.suffix](path, dtype=str) for gene_type, path in cc_sr_paths.items()
    }

    # The two currently supported species have different file types and
    # formats, so they need to be handled differently
    if 'hsapiens' in ds:
        # For humans, the two gene types are in CSVs that are read into
        # pandas dataframes. Get a series of gene IDs for each
        gene_ids_by_type = {
            'cell_cycle': cc_sr_genes['cell_cycle'].loc[
                cc_sr_genes['cell_cycle']['Manual_annotation'] != 'Uncharacterised', 'Ensembl_ID'
            ],
            'stress_response': cc_sr_genes['stress_response'].loc[:, 'ensembl_gene_id'],
        }

        # Iterate over the two gene types and add as boolean columns to
        # the genes_df
        for gene_type, gene_ids in gene_ids_by_type.items():
            genes_df[gene_type] = genes_df['ensembl_gene_id'].isin(gene_ids)

    elif 'mmusculus' in ds:
        # For mice, they are just two plain text files containing gene
        # names, so assign those
        for gene_type, gene_ids in cc_sr_genes.items():
            genes_df[gene_type] = genes_df['external_gene_name'].isin(gene_ids)

    # Get a dataframe of the datasets in this mart
    dsets_df = s.get_datasets(mart)

    # Get the description of the one dataset matching the current dataset
    description = dsets_df.loc[dsets_df['name'] == ds, 'description'].values[0]

    # The genome is encapsulated in parentheses inside of this description
    # so use re.search to get it out
    genome = search(pattern=r'\((.*)\)', string=description).group(1)

    # Finally, save to parquet file for easy reading in pandas and R
    output_path = Path(f'../assets/ref_annotations/{genome}.parquet')
    genes_df.to_parquet(output_path)
