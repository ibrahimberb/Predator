from collections import defaultdict
import pandas as pd


def get_patient_protein_to_mutations_dict(patient_snv_data: pd.DataFrame):
    """  
    Given a patient_snv_data frame, which is filtered for a specific patient ID, return
    the dictionary `patient_protein_to_mutations` that maps protein to its mutations.
    E.g.
        'Q8TDI0': ['L33V'],
        'P15812': ['I163F', 'I163N']
        ...

    Parameters
    ----------
        patient_snv_data : <DataFrame>
            A SNV data which contains protein and mutations for only one patient.

    Returns
    -------
        patient_protein_to_mutations : <defaultdict>
            A (default) dictionary that maps each protein to its mutation positions.

    """
    
    # Confirm that given patient_snv_data is contains only one patient id.
    assert len(patient_snv_data['Tumor_Sample_Barcode'].unique()) == 1

    # Dictionary that maps protein <str> to mutations <list> for current patient.
    patient_protein_to_mutations = defaultdict(list)

    # Fill in the dictionary `patient_protein_to_mutations`.
    for index, row in patient_snv_data.iterrows():
        protein = row['SWISSPROT']
        mutation = row['HGVSp_Short']
        patient_protein_to_mutations[protein].append(mutation)

    return patient_protein_to_mutations
