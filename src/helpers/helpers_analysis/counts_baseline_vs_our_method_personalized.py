import pandas as pd
from tqdm.notebook import tqdm

from ..labels import ClassLabels

from .get_patient_protein_to_mutations_dict import get_patient_protein_to_mutations_dict
from .is_core import is_core
from .is_in_elaspic import is_in_elaspic

from collections import defaultdict


def counts_baseline_vs_our_method_personalized(
        proteins: list, patients: list,
        snv_data: pd.DataFrame,
        elaspic_core_data: pd.DataFrame,
        elaspic_interface_data: pd.DataFrame, prediction_data: pd.DataFrame,
        add_core_flag_1_case_dict: dict
):
    """
    Generates personalized protein counts for BASELINE and OUR_METHOD.

    Parameters
    ----------
        proteins : <list>
            List of proteins.

        patients : <list>
            A list of TCGA patients.

        snv_data : <DataFrame>
            An SNV dataframe, we use the processed version of SNV.

        elaspic_core_data : <DataFrame>
            The ELASPIC results file that contains only the `core` type entries.

        elaspic_interface_data : <DataFrame>
            The ELASPIC results file that contains only the `interface` type entries. It will be used to
            check if a specific (protein, mutation) pair is an interface via `is_interface` function.

        prediction_data : <DataFrame>
            The dataframe which contains prediction column, along with protein, mutation, interactor columns.
        
        add_core_flag_1_case_dict : <None> or <dict>
            Controls whether to add `ELASPIC degree` or to add +0 (i.e. ignoring `core_flag=1` case).
            If `None`, `core_flag=1` case adds +0.
            Otherwise `core_flag=1` case adds ELASPIC degree.

    Returns
    -------
        proteins_to_counts_baseline_dict : <dict>
            A dictionary that maps each protein to counts for BASELINE.

        proteins_to_counts_our_method_dict : <dict>
            A dictionary that maps each protein to counts for OUR_METHOD.

    """

    if add_core_flag_1_case_dict:
        print('Adding ELASPIC Degree when `core_flag=1`')  # TODO

    else:
        print('Adding +0 when `core_flag=1`')
        add_core_flag_1_case_dict = defaultdict(int)

    # MAIN: PERSONALIZED
    personalized_proteins_to_counts_baseline_dict = dict()
    personalized_proteins_to_counts_our_method_dict = dict()

    for patient in tqdm(patients):
        # Setting the counts 0.
        proteins_to_counts_baseline_dict = dict.fromkeys(proteins, 0)
        proteins_to_counts_our_method_dict = dict.fromkeys(proteins, 0)

        # Filter SNV data for current patient.
        patient_snv_data = snv_data[snv_data["Tumor_Sample_Barcode"] == patient]

        for protein, mutations in get_patient_protein_to_mutations_dict(patient_snv_data).items():

            core_flag = 'N/A'
            # print(protein, mutations)
            for mutation in mutations:

                # Check if (protein.mutation) is in ELASPIC.
                if is_in_elaspic(protein, mutation, elaspic_core_data, elaspic_interface_data):
                    # print(f'{protein}.{mutation} IS IN ELASPIC.')

                    if is_core(protein, mutation, elaspic_core_data):
                        # print(' → core found!')
                        core_flag = 1
                        break

                    else:
                        # print(' → interface found!')
                        core_flag = 0

                else:
                    # print(f'{protein}.{mutation} IS NOT IN ELASPIC.')
                    # print(f'CORE_FLAG = {core_flag}')
                    continue

            if core_flag == 1:
                # print(f'CORE_FLAG = {core_flag}')
                # print('+ Adding counts.. ', brca_proteins_to_elaspic_degree[protein])
                # increase baseline counts by elaspic degree
                # Increase our method counts by elaspic degree
                proteins_to_counts_baseline_dict[protein] += add_core_flag_1_case_dict[protein]
                proteins_to_counts_our_method_dict[protein] += add_core_flag_1_case_dict[protein]

            elif core_flag == 0:
                # For our model: Contains Disruptive interactions only.
                disruptive_interactions = set()
                for mutation in mutations:
                    # Check if (protein.mutation) is in ELASPIC.
                    if is_in_elaspic(protein, mutation, elaspic_core_data, elaspic_interface_data):
                        prediction_search = prediction_data[
                            (prediction_data['UniProt_ID'] == protein) &
                            (prediction_data['Mutation'] == mutation) &
                            (prediction_data['Prediction'] == ClassLabels.DISRUPTING)].copy()
                        # add interactor proteins to disruptive_interactions set.
                        interactor_list = prediction_search['Interactor_UniProt_ID'].to_list()
                        disruptive_interactions.update(interactor_list)

                # For baseline model: Contains Increasing+NoEff interactions and Disruptive interactions
                interactions = set()
                for mutation in mutations:
                    # Check if (protein.mutation) is in ELASPIC.
                    if is_in_elaspic(protein, mutation, elaspic_core_data, elaspic_interface_data):
                        prediction_search = prediction_data[
                            (prediction_data['UniProt_ID'] == protein) &
                            (prediction_data['Mutation'] == mutation)].copy()
                        # add interactor proteins to interactions set.
                        interactor_list = prediction_search['Interactor_UniProt_ID'].to_list()
                        interactions.update(interactor_list)

                # print(f'CORE_FLAG = {core_flag}')
                # print('+ Adding counts.. ')
                # increase baseline counts by elaspic degree
                proteins_to_counts_baseline_dict[protein] += len(interactions)
                # Increase our method counts by depending our predictions
                proteins_to_counts_our_method_dict[protein] += len(disruptive_interactions)

        personalized_proteins_to_counts_baseline_dict[patient] = proteins_to_counts_baseline_dict
        personalized_proteins_to_counts_our_method_dict[patient] = proteins_to_counts_our_method_dict

    return personalized_proteins_to_counts_baseline_dict, personalized_proteins_to_counts_our_method_dict
