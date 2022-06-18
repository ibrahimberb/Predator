from tqdm.notebook import tqdm
from collections import defaultdict
import pandas as pd

from .get_patient_protein_to_mutations_dict import get_patient_protein_to_mutations_dict
from .is_core import is_core
from .is_in_elaspic import is_in_elaspic


def get_patients_to_interface_proteins_dict(
        patients: list, snv_data: pd.DataFrame,
        elaspic_core_data: pd.DataFrame,
        elaspic_interface_data: pd.DataFrame,
        print_details=False
):
    """
    Given a patient list, return the dictionary `patients_to_interface_proteins` that maps each patient
    to its interface only proteins, i.e. proteins whose `core_flag = 0`.

    Parameters
    ----------
        # TODO: DOC Strings..
        snv_data : <DataFrame>
            An SNV dataframe, we use the processed version of SNV.

        elaspic_core_data : <DataFrame>
            The ELASPIC results file that contains only the `core` type entries.

        elaspic_interface_data : <DataFrame>
            The ELASPIC results file that contains only the `interface` type entries.  

    Returns
    -------
        patients_to_interface_proteins : <defaultdict>
            A (default) dictionary that maps each patient to its interface_proteins.

    """

    # Dictionary that maps patients <str> to interface_proteins <list> for current patient.
    patients_to_interface_proteins = defaultdict(list)

    for patient in tqdm(patients):
        if print_details: print(f'Current patient: {patient}')
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
                pass

            elif core_flag == 0:
                if print_details: print('PROTEIN:', protein)
                # adding to the `brca_patients_to_interface_proteins` dictionary.
                patients_to_interface_proteins[patient].append(protein)
                # print(f'CORE_FLAG = {core_flag}')
                # print('+ Adding counts.. ')

    return patients_to_interface_proteins
