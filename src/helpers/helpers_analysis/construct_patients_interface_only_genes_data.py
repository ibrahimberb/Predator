import pandas as pd


def construct_patients_interface_only_genes_data(
    patients: list, 
    patients_to_interface_proteins: dict,
    protein_to_gene_dict: dict, 
    cgc_type_to_genes: dict
):
    """
    Constructs a dataframe that stores interface only genes (i.e. `core_flag=0`)
    for each patient, and their intersection with CGC genes.

    Parameters
    ----------
        patients : <list>
            A list of TCGA patients.

        patients_to_interface_proteins : <dict>
            A (default) dictionary that maps each patient to its interface_proteins.

        protein_to_gene_dict : <dict>
            A dictionary that maps each protein to corresponding gene, used for converting protein
            names into gene names.

        cgc_type_to_genes : <dict>
            A dictionary that holds list of CGC genes based on their cohort type.

    Returns
    -------
        patients_to_interface_proteins : <defaultdict>
            A (default) dictionary that maps each patient to its interface_proteins.
    """

    # Create the dataframe `patients_interface_only_genes_data` using `patients`.
    patients_interface_only_genes_data = pd.DataFrame(index=patients)

    # Adding the `Interface_only_proteins` column.
    interface_only_proteins = []
    for patient in patients:
        interface_only_proteins.append(patients_to_interface_proteins[patient])

    patients_interface_only_genes_data['Interface_only_proteins'] = interface_only_proteins

    # Adding the `Interface_only_genes` column: Converting proteins in `Interface_only_proteins` to genes.
    interface_only_genes = []
    for patient in patients:
        # For current patient, convert protein list to their corresponding genes. E.g. ['P1', 'P2'] → ['G1', 'G2']
        interface_only_genes.append([protein_to_gene_dict[protein] for protein in patients_to_interface_proteins[patient]])

    patients_interface_only_genes_data['Interface_only_genes'] = interface_only_genes

    # Adding number of proteins.
    patients_interface_only_genes_data['Num_Interface_only_genes'] = patients_interface_only_genes_data['Interface_only_genes'].apply(len)
    
    # Adding CGC Intersection
    cgc_intersection = []
    for patient in patients:
        current_interface_only_genes = [protein_to_gene_dict[protein] for protein in patients_to_interface_proteins[patient]]
        cgc_intersection.append([gene for gene in current_interface_only_genes if gene in cgc_type_to_genes['CGC']])

    patients_interface_only_genes_data['CGC_intersection'] = cgc_intersection

    # Adding Num_cgc_intersection
    patients_interface_only_genes_data['Num_CGC_intersection'] = patients_interface_only_genes_data['CGC_intersection'].apply(len)
    
    # Adding CGC BRCA Intersection
    cgc_brca_intersection = []
    for patient in patients:
        current_interface_only_genes = [protein_to_gene_dict[protein] for protein in patients_to_interface_proteins[patient]]
        cgc_brca_intersection.append([gene for gene in current_interface_only_genes if gene in cgc_type_to_genes['CGC_BRCA']])

    patients_interface_only_genes_data['CGC_brca_intersection'] = cgc_brca_intersection

    # Adding Num_CGC_intersection
    patients_interface_only_genes_data['Num_CGC_brca_intersection'] = patients_interface_only_genes_data['CGC_brca_intersection'].apply(len)

    # Return the constructed dataframe.
    return patients_interface_only_genes_data
