import pandas as pd


def get_single_patient_personalized_data(
        patients_interface_only_genes_data: pd.DataFrame, protein_to_gene_dict: dict,
        personalized_proteins_to_counts_baseline_param: dict,
        personalized_proteins_to_counts_our_method_param: dict,
        patient_id: str
):
    """

    Parameters
    ----------
    patients_interface_only_genes_data
    protein_to_gene_dict
    personalized_proteins_to_counts_baseline_param
    personalized_proteins_to_counts_our_method_param
    patient_id

    Returns
    -------
    personalized_data
    """

    personalized_data = pd.DataFrame(patients_interface_only_genes_data.loc[patient_id, 'Interface_only_genes'],
                                     columns=['GENE'])

    personalized_baseline = []
    personalized_our_method = []
    for gene in patients_interface_only_genes_data.loc[patient_id, 'Interface_only_genes']:
        # converting gene name to protein name.
        converted_protein = list(protein_to_gene_dict.keys())[list(protein_to_gene_dict.values()).index(gene)]
        personalized_baseline.append(personalized_proteins_to_counts_baseline_param[patient_id][converted_protein])
        personalized_our_method.append(personalized_proteins_to_counts_our_method_param[patient_id][converted_protein])

    personalized_data['BASELINE'] = personalized_baseline
    personalized_data['OUR_METHOD'] = personalized_our_method
    personalized_data['OUR_METHOD_NORMALIZED'] = personalized_data['OUR_METHOD'] / personalized_data['BASELINE']

    return personalized_data
