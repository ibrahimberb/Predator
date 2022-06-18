from tqdm.notebook import tqdm
import pandas as pd

from .get_patient_protein_to_mutations_dict import get_patient_protein_to_mutations_dict
from ..labels import ClassLabels
from .is_core import is_core
from .is_in_elaspic import is_in_elaspic


def add_baseline(preliminary_data, baseline_dict):
    """
    Given a tcga preliminary_data data, add baseline counts.

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `BASELINE` column will be to be added on.

        baseline_dict : <dict>
            A dictionary that maps proteins to baseline counts.

    Returns
    -------
        None. Modifies the dataframe.
    """

    baseline_counts = []
    for protein in preliminary_data["PROTEIN"]:
        baseline_counts.append(baseline_dict[protein])

    preliminary_data['BASELINE'] = baseline_counts


def add_cancermine_status(preliminary_data, cancermine_genes, cancermine_cohort_genes, tcga_type):
    """
    Given a tcga preliminary_data data, add the CancerMine and CancerMine_Cohort status. In other words, find
    if given gene is CancerMine and CancerMine_Cohort gene lists.
    E.g.
        ----+------------+-------------------+
        ..  | CancerMine | CancerMine_Cohort |
        ----+------------+-------------------+
        ..  |    YES     |         NO        |
        ..  |     NO     |         NO        |
        ..  |     ..     |         ..        |
        ----+------------+-------------------+


    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `CancerMine_STATUS` and `CancerMine_Cohort_STATUS` column will be to be added on.

        cancermine_genes : <list>
            CancerMine gene list.

        cancermine_cohort_genes : <list>
            CancerMine cohort gene list. Must be one of BRCA, COAD, or OV, depending on `preliminary_data`.

        tcga_type : <str>
            The TCGA abbreviation, used for column name. Must be one of BRCA, COAD, or OV.

    Returns
    -------
        None. Modifies the dataframe.
    """
    cancermine_statuses = []
    cancermine_cohort_statuses = []

    for gene in preliminary_data["GENE"]:
        cancermine_status = "+" if gene in cancermine_genes else "-"
        cancermine_cohort_status = "+" if gene in cancermine_cohort_genes else "-"
        cancermine_statuses.append(cancermine_status)
        cancermine_cohort_statuses.append(cancermine_cohort_status)

    # Add the columns
    preliminary_data['CancerMine_STATUS'] = cancermine_statuses
    preliminary_data[f'CancerMine_STATUS ({tcga_type})'] = cancermine_cohort_statuses


def add_cgc_status(preliminary_data, cgc_genes, cgc_cohort_genes, tcga_type):
    """
    Given a tcga preliminary_data data, add the CGC and CGC_Cohort status. In other words, find 
    if given gene is CGC and CGC_Cohort gene lists.
    E.g.
        ----+-----+------------+
        ..  | CGC | CGC_Cohort |
        ----+-----+------------+
        ..  | YES |     NO     |
        ..  | NO  |     NO     |
        ..  | ..  |     ..     |
        ----+-----+------------+


    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `CGC_STATUS` and `CGC_Cohort_STATUS` column will be to be added on.

        cgc_genes : <list>
            CGC gene list.

        cgc_cohort_genes : <list>
            CGC cohort gene list. Must be one of BRCA, COAD, or OV, depending on `preliminary_data`.

        tcga_type : <str>
            The TCGA abbreviation, used for column name. Must be one of BRCA, COAD, or OV.

    Returns
    -------
        None. Modifies the dataframe.
    """
    cgc_statuses = []
    cgc_cohort_statuses = []

    for gene in preliminary_data["GENE"]:
        cgc_status = "+" if gene in cgc_genes else "-"
        cgc_cohort_status = "+" if gene in cgc_cohort_genes else "-"
        cgc_statuses.append(cgc_status)
        cgc_cohort_statuses.append(cgc_cohort_status)

    # Add the columns
    preliminary_data['CGC_STATUS'] = cgc_statuses
    preliminary_data[f'CGC_STATUS ({tcga_type})'] = cgc_cohort_statuses


def add_elaspic_coverage(preliminary_data, elaspic_combined_data, snv_data):
    """
    Given preliminary_data, elaspic_combined_data and snv_data;
    
        Go over each patients in SNV data.
        List all specific protein.mutation pairs for current patient.
        If *any* of these mutations are in ELASPIC, increase the count (without considering CORE or INTERFACE)
        Then, go to the next patient.
        
    This can be thought of as another kind of baseline.
        
    Note: We count only once if there are multiple occurrance of the same protein's mutations.
    
    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `ELASPIC_COVERAGE` column will be to be added on.
            
        elaspic_combined_data : <dataframe>
            The ELASPIC combined results file which contains both CORE and INTERFACE.
        
        snv_data : <DataFrame>
            An SNV dataframe, we use the processed version of SNV.
    
    Returns
    -------
        None. Modifies the dataframe.
    """
    
    # Initialize the `proteins_to_counts_dict` dictionary with proteins in preliminary_data in corrent order, 
    # and the default values of 0.
    proteins_to_counts_dict = dict.fromkeys(list(preliminary_data['PROTEIN']), 0)
    
    # Get the patient IDs.
    patients = list(snv_data['Tumor_Sample_Barcode'].unique())
    
    for patient in tqdm(patients):
        
        # Patient filtered dataframe: Filter SNV file for current patient.
        patient_snv_data = snv_data[snv_data["Tumor_Sample_Barcode"] == patient]
        
        # Get the dictionary that maps protein <str> to mutations <list> for *current patient*.
        patient_protein_to_mutations = get_patient_protein_to_mutations_dict(patient_snv_data)      
        
        # Fill the counts dictionary `proteins_to_counts_dict`
        for protein, mutations in patient_protein_to_mutations.items():
            for mutation in mutations:
                
                # Search in ELASPIC.
                elaspic_search_data = elaspic_combined_data[(elaspic_combined_data["UniProt_ID"] == protein) & 
                                                            (elaspic_combined_data["Mutation"] == mutation)]
                
                # If search data is not empty, i.e. A specific Mutation of current Protein
                # exits in ELASPIC data, increase the count and skip to the next Protein.
                if not elaspic_search_data.empty:
                    proteins_to_counts_dict[protein] += 1
                    break

                    
    # Add the column
    preliminary_data['ELASPIC_COVERAGE'] = list(proteins_to_counts_dict.values())


def add_genes(preliminary_data, protein_to_gene_dict):
    """
    Given a tcga preliminary_data data, add corresponding genes.
    E.g.
        P62837 → UBE2D2

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `GENES` column will be to be added on.

        protein_to_gene_dict : <dict>
            A dictionary that mappes protein name to gene name. 
            One of the followings:
            - brca_protein_to_gene_dict
            - coad_protein_to_gene_dict
            - ov_protein_to_gene_dict

    Returns
    -------
        None. Modifies the dataframe.
    """

    genes = []
    for protein in preliminary_data["PROTEIN"]:
        genes.append(protein_to_gene_dict[protein])

    preliminary_data['GENE'] = genes


def add_num_disruptive_entries(preliminary_data: pd.DataFrame, prediction_data: pd.DataFrame):
    """
    Given a tcga preliminary_data data, find the number of disruptive predicted entries for each protein.
    E.g.
        P04637 (TP53) → 76 → 74
        
    Note: For some entries, it is possible that we may not have any predictions, either because it is 
        not an interface mutation or it was causing "invalid" predictions with its triplet 
        (Disrupting and nondisrupting at the same time).

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `NUM_DISRUPTIVE_ENTRIES` column will be to be added on.

        prediction_data : <DataFrame>
            The dataframe which contains prediction column, along with protein, mutation, interactor columns.

    Returns
    -------
        None. Modifies the dataframe.
    """

    num_disruptive_entries = []
    for protein in preliminary_data["PROTEIN"]:
        # Get the disruptive entries
        disruptive_entries_data = prediction_data[(prediction_data['UniProt_ID'] == protein) &
                                                  (prediction_data['Prediction'] == ClassLabels.DISRUPTING)]
        num_disruptive_entries.append(len(disruptive_entries_data))

    preliminary_data['NUM_DISRUPTIVE_ENTRIES'] = num_disruptive_entries


def add_num_incr_noeff_entries(preliminary_data: pd.DataFrame, prediction_data: pd.DataFrame):
    """
    Given a tcga preliminary_data data, find the number of increasing+no_effect predicted entries for each protein.
    E.g.
        P04637 (TP53) → 76 → 1
        
    Note: For some entries, it is possible that we may not have any predictions, either because it is 
        not an interface mutation or it was causing "invalid" predictions with its triplet
        (Disrupting and nondisrupting at the same time).

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `NUM_INCR_NOEFF_ENTRIES` column will be to be added on.

        prediction_data : <DataFrame>
            The dataframe which contains prediction column, along with protein, mutation, interactor columns.

    Returns
    -------
        None. Modifies the dataframe.
    """

    num_incr_noeff_entries = []
    for protein in preliminary_data["PROTEIN"]:
        # Get the increasing + no_effect entries (nondisrupting).
        incr_noeff_entries_data = prediction_data[(prediction_data['UniProt_ID'] == protein) &
                                                  (prediction_data['Prediction'] == ClassLabels.NONDISRUPTING)]
        num_incr_noeff_entries.append(len(incr_noeff_entries_data))

    preliminary_data['NUM_INCR_NOEFF_ENTRIES'] = num_incr_noeff_entries


def add_num_interface_patients_disruptive_interactor(
        preliminary_data, patient_interaction_data
):
    """
    Given a tcga preliminary_data and patient interaction data for that cohort,
    add number of interface patients where *at least one* interaction is disrupted
    for that gene.

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `NUM_INTERFACE_PATIENTS_DISRUPTIVE_INTERACTOR` column
            will be added on.

        patient_interaction_data : <DataFrame>
            The patient interaction data with disruptive interactors for each patient. Also, there is
            CORE / INTERFACE categorization for each patient, shown with "C" or "I".

    Returns
    -------
        None. Modifies the dataframe.
    """

    num_interface_patients_disruptive_interactor = []
    for protein in preliminary_data["PROTEIN"]:
        num_patients = patient_interaction_data[
            # Looking for Interface patients only -- with "I".
            (patient_interaction_data["CORE_INTERFACE_VS_INTERFACE_STATUS"] == "I") &
            # Looking for at least one disruptive interactors.
            (patient_interaction_data["NUM_DISRUPTIVE_INTERACTORS"] > 0) &
            # Looking for current protein.
            (patient_interaction_data["PROTEIN_GENE"].apply(lambda x: x.split(":")[0] == protein))
            ]["PATIENT"].nunique()

        num_interface_patients_disruptive_interactor.append(num_patients)

    preliminary_data['NUM_INTERFACE_PATIENTS_DISRUPTIVE_INTERACTOR'] = num_interface_patients_disruptive_interactor


def add_num_elaspic_interface_entries(preliminary_data, protein_to_num_elaspic_interface_entries):
    """
    Given a tcga preliminary_data data, add corresponding number of ELASPIC interface entries for each protein.
    E.g.
        P04637 (TP53) → 76

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `NUM_ELASPIC_INTERFACE_ENTRIES` column will be to be added on.

        protein_to_num_elaspic_interface_entries : <defaultdict>
            A (default) dictionary that maps each protein to its number of ELASPIC interface entries.

    Returns
    -------
        None. Modifies the dataframe.
    """

    elaspic_interface_entries = []
    for protein in preliminary_data["PROTEIN"]:
        elaspic_interface_entries.append(protein_to_num_elaspic_interface_entries[protein])

    preliminary_data['NUM_ELASPIC_INTERFACE_ENTRIES'] = elaspic_interface_entries


def add_num_unique_interactors(preliminary_data, protein_to_num_unique_interactors):
    """                                          
    Given a tcga preliminary_data data, add corresponding number of unique ELASPIC interface interactors 
    for each protein.
    E.g.
        P04637 (TP53) → 15

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `NUM_UNIQUE_INTERACTORS` column will be to be added on.

        protein_to_num_unique_interactors : <defaultdict>
            A (default) dictionary that maps each protein to its number of unique interactors in ELASPIC interface data.


    Returns
    -------
        None. Modifies the dataframe.
    """

    num_unique_interactors = []
    for protein in preliminary_data["PROTEIN"]:
        num_unique_interactors.append(protein_to_num_unique_interactors[protein])

    preliminary_data['NUM_UNIQUE_INTERACTORS'] = num_unique_interactors


def add_our_method(preliminary_data, our_method_dict):
    """
    Given a tcga preliminary_data data, add our_method counts.

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `OUR_METHOD` column will be to be added on.

        our_method_dict : <dict>
            A dictionary that maps proteins to our_method counts.

    Returns
    -------
        None. Modifies the dataframe.
    """

    our_method_counts = []
    for protein in preliminary_data["PROTEIN"]:
        our_method_counts.append(our_method_dict[protein])

    preliminary_data['OUR_METHOD'] = our_method_counts


def add_patient_core_count(preliminary_data, snv_data, elaspic_core_data, elaspic_interface_data):
    """
    Given a tcga preliminary_data data, add number of patients that core mutation occurred
    for each protein as a new column.

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `PATIENT_CORE_COUNT` column will be to be added on.

        snv_data : <DataFrame>
            An SNV dataframe, we use the processed version of SNV.

        elaspic_core_data : <DataFrame>
            The ELASPIC results file that contains only the `core` type entries.

        elaspic_interface_data : <DataFrame>
            The ELASPIC results file that contains only the `interface` type entries. It will be used to
            check if a specific (protein, mutation) pair is an interface via `is_interface` function.

    Returns
    -------
        None. Modifies the dataframe.
    """

    # Initialize the `proteins_to_patient_core_counts_dict` dictionary with proteins in preliminary_data in
    # correct order, and the default values of 0.
    proteins_to_patient_core_counts_dict = dict.fromkeys(list(preliminary_data['PROTEIN']), 0)

    # Get the patient IDs.
    patients = list(snv_data['Tumor_Sample_Barcode'].unique())

    for patient in tqdm(patients):

        # Patient filtered dataframe: Filter SNV file for current patient.
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
                # Adding the corresponding gene counter
                proteins_to_patient_core_counts_dict[protein] += 1

    # Add the column
    preliminary_data['PATIENT_CORE_COUNT'] = list(proteins_to_patient_core_counts_dict.values())


def add_patient_interface_count(preliminary_data, snv_data, elaspic_core_data, elaspic_interface_data):
    """
    Given a tcga preliminary_data data, add number of patients that only interface mutation occurred
    for each protein as a new column.

    Parameters
    ----------
        preliminary_data : <DataFrame>
            Preliminary data which `PATIENT_INTERFACE_COUNT` column will be to be added on.

        snv_data : <DataFrame>
            An SNV dataframe, we use the processed version of SNV.

        elaspic_core_data : <DataFrame>
            The ELASPIC results file that contains only the `core` type entries.

        elaspic_interface_data : <DataFrame>
            The ELASPIC results file that contains only the `interface` type entries. It will be used to
            check if a specific (protein, mutation) pair is an interface via `is_interface` function.

    Returns
    -------
        None. Modifies the dataframe.
    """

    # Initialize the `proteins_to_patient_interface_counts_dict` dictionary with proteins in preliminary_data in
    # correct order, and the default values of 0.
    proteins_to_patient_interface_counts_dict = dict.fromkeys(list(preliminary_data['PROTEIN']), 0)

    # Get the patient IDs.
    patients = list(snv_data['Tumor_Sample_Barcode'].unique())

    for patient in tqdm(patients):

        # Patient filtered dataframe: Filter SNV file for current patient.
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

            if core_flag == 0:
                # Adding the corresponding gene counter
                proteins_to_patient_interface_counts_dict[protein] += 1

    # Add the column
    preliminary_data['PATIENT_INTERFACE_COUNT'] = list(proteins_to_patient_interface_counts_dict.values())

