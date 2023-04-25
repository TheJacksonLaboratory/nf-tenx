import pandas as pd
from bioservices import BioMart
from io import StringIO
from pathlib import Path
from numpy import loadtxt

ensembl_datasets = ["hsapiens_gene_ensembl", "mmusculus_gene_ensembl"]
attributes_of_interest = [
    "ensembl_gene_id",
    "external_gene_name",
    "chromosome_name",
    "description",
]

for ds in ensembl_datasets:
    # Establish a connection to ENSEMBL biomart and add the dataset to
    # the XML query
    s = BioMart()
    s.add_dataset_to_xml(ds)

    # Add each attribute to the XML query
    for attribute in attributes_of_interest:
        s.add_attribute_to_xml(attribute)

    # Generate the XML query and store the result, convert to buffer
    xml_query = s.get_xml()
    result_string = s.query(xml_query)
    result_buffer = StringIO(result_string)

    # Use buffer to make a dataframe, where genes will be annotated
    genes_df = pd.read_csv(result_buffer, sep="\t", names=attributes_of_interest)

    # Set a boolean column for each gene type using some condition

    # Mitochondrial genes are on the mitochondrial chromosome
    genes_df["mitochondrial"] = genes_df["chromosome_name"].str.match("mt", case=False)

    # Assume that if hemoglobin is in the description, then it's a
    # hemoglobin gene. Also, filter out NaNs
    genes_df["hemoglobin"] = (genes_df["description"].str.contains("hemoglobin", case=False)) & (genes_df["description"] == genes_df["description"])

    # Sex linked genes are either on the Y-chromosome or are the
    # x-inactivation gene
    genes_df["sex_linked"] = (genes_df["chromosome_name"].str.match("y", case=False)) | (genes_df["external_gene_name"].str.contains("xist", case=False))

    # Ribosomal gene names have a certain pattern, so match that
    # pattern. Also filter out NaNs because not all genes have external
    # gene names
    genes_df["ribosomal"] = (genes_df["external_gene_name"].str.match("m?rp[ls]", case=False)) & (genes_df["external_gene_name"] == genes_df["external_gene_name"])

    # Define directory for species, where gene CSVs live
    species_dir = Path(ds)

    # The two currently supported species have different file types
    # and formats, so they need to be handled differently
    if ds == "hsapiens_gene_ensembl":
        # Read cell cycle genes into pandas dataframe
        cell_cycle_path = species_dir / Path("cell_cycle_genes.xlsx")
        cell_cycle_df = pd.read_excel(cell_cycle_path)

        # Get a series of gene_ids that for characterized genes
        cell_cycle_genes = cell_cycle_df.loc[cell_cycle_df["Manual_annotation"] != "Uncharacterised", "Ensembl_ID"]

        # Do the same for stress response genes. Note that this CSV has no
        # "Manual Annotation" column, so just take all of them
        stress_response_path = species_dir / Path("stress_response_genes.csv")  
        stress_response_df = pd.read_csv(stress_response_path)
        stress_response_genes = stress_response_df.loc[:, "ensembl_gene_id"]

        # Finally, assign a column of the genes_df for all gene_ids in
        # these series
        genes_df["cell_cycle"] = genes_df["ensembl_gene_id"].isin(cell_cycle_genes)
        genes_df["stress_response"] = genes_df["ensembl_gene_id"].isin(stress_response_genes)

    elif ds == "mmusculus_gene_ensembl":
        # Mouse cell cycle gene names are in a plain text file, use nup
        # to load them in
        cell_cycle_path = species_dir / Path("cell_cycle_genes.txt")
        cell_cycle_genes = loadtxt(cell_cycle_path, dtype=str)

        # Same for stress response genes
        stress_response_path = species_dir / Path("stress_response_genes.txt")
        stress_response_genes = loadtxt(stress_response_path, dtype=str)
        
        # Assign column of genes_df according to gene name instead of
        # ensembl_gene_id
        genes_df["cell_cycle"] = genes_df["external_gene_name"].isin(cell_cycle_genes)
        genes_df["stress_response"] = genes_df["external_gene_name"].isin(stress_response_genes)

    # Finally, save to pickled dataframe object
    genes_df.to_pickle(species_dir / Path("annotated_genes"))