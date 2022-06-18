# Disruptive Interactions for Each Patients.

# read the prediction data

# read the SNV data and get `patient_to_mutations` dictionary.

# For each patient, obtain the disruptive interactions located in prediction.
# We look for protein.mutation only here.

from collections import defaultdict
from datetime import datetime
from typing import List

from pandas import read_csv
import pandas as pd
from pathlib import Path
import os.path as op

from tqdm.auto import tqdm
# import os.path as op

from ..mylogger import get_handler
import logging

from .gene_id_retrieval import GeneIDRetriever

from .get_patient_protein_to_mutations_dict import get_patient_protein_to_mutations_dict

from .loaders import (
    load_snv_datasets,
    load_elaspic_datasets
)

from ..labels import ClassLabels

from .is_core import is_core
from .is_in_elaspic import is_in_elaspic

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class InteractionsPerPatient:
    def __init__(
            self,
            tcga: str,
            prediction_data_path: Path,
            tcga_snv_path: Path,
            elaspic_core_path: Path,
            elaspic_interface_path: Path,
            identifier: str,
            verbose: bool = False
    ):
        """
        Interactions per Patient.

        Parameters
        ----------
            tcga : <str>
                TCGA Cohort name, e.g. `BRCA`.

            prediction_data_path : <Path>
                The path of prediction data.

            tcga_snv_path : <Path>
                The path of original TCGA SNV data.

            identifier : <str>
                The identifier to be used in protein/gene representation.
                - `uniprot` for UniProt_ID protein representation.
                - `hugo` for Hugo ID Gene representation.

            verbose : <bool>
                Controls the detail progress to be displayed. E.g. decides whether patient in
                progress to be displayed when finding that patient's disruptive interactions.

        """
        self.tcga = tcga.upper()
        self.prediction_data_path = prediction_data_path
        self.tcga_snv_path = tcga_snv_path
        self.elaspic_core_path = elaspic_core_path
        self.elaspic_interface_path = elaspic_interface_path
        self.identifier = identifier
        self.verbose = verbose

        self.snv_data_simplified = None
        self.patients = None
        self.patient_to_snv_data = None
        self.disruptive_prediction_data = None
        self.prediction_data = None
        self.analysis_table = None

        self.core_data = None
        self.interface_data = None

        self.patient_to_interactions = None
        self.patient_to_disruptive_interactions = None

        self.patient_protein_to_C_I_status_data = None

        self.uniprot_to_gene_id = None
        self.gene_retriever = GeneIDRetriever()
        self.load_materials()
        self.prepare_patient_protein_to_C_I_status_data(self.patients)

        self.find_all_interactions()
        self.find_disruptive_interactions()

        log.info("Disruptive interactions per patient completed.")

    def load_materials(self):
        log.info("Loading materials ..")
        self.load_snv_data_simplified()
        self.load_patient_ids()
        self.load_patient_to_snv_data()
        self.load_prediction_data()
        self.load_disruptive_prediction_data()
        self.load_uniprot_to_gene_id()
        self.load_elaspic_core_and_interface_datasets()
        log.info("Materials loaded.")
        log.info(f"Number of {self.tcga} patients: {len(self.patients)}.")

    def load_elaspic_core_and_interface_datasets(self):
        elaspic_data_materials = {}
        load_elaspic_datasets(
            tcga=self.tcga,
            elaspic_core_path=self.elaspic_core_path,
            elaspic_interface_path=self.elaspic_interface_path,
            data_materials=elaspic_data_materials,
        )

        self.core_data = elaspic_data_materials[f"{self.tcga}_elaspic_core_data_simplified"]
        self.interface_data = elaspic_data_materials[f"{self.tcga}_elaspic_interface_processed_data"]

    def load_snv_data_simplified(self):
        log.debug("Loading SNV data simplified ..")
        data_materials = {}
        load_snv_datasets(self.tcga, self.tcga_snv_path, data_materials)
        self.snv_data_simplified = data_materials[f"{self.tcga}_snv_data_simplified"]

    def load_patient_ids(self):
        log.debug("Loading patient ids ..")
        patients = list(self.snv_data_simplified["Tumor_Sample_Barcode"].unique())
        self.patients = patients

    def get_patient_snv_data(self, patient_id):
        patient_snv_data = self.snv_data_simplified[
            self.snv_data_simplified["Tumor_Sample_Barcode"] == patient_id
            ]
        return patient_snv_data

    def load_patient_to_snv_data(self):
        log.debug("Loading patient to snv_data ..")
        patient_to_snv_data = {}
        for patient in tqdm(self.patients):
            patient_snv_data = self.get_patient_snv_data(patient)
            patient_to_snv_data[patient] = patient_snv_data

        self.patient_to_snv_data = patient_to_snv_data

    def load_prediction_data(self):
        log.debug("Loading the prediction data ..")
        prediction_data = read_csv(self.prediction_data_path)
        self.prediction_data = prediction_data

    def load_disruptive_prediction_data(self):
        log.debug("Loading the disruptive prediction data ..")
        disruptive_prediction_data = self.prediction_data[
            self.prediction_data["Prediction"] == ClassLabels.DISRUPTING
        ]
        self.disruptive_prediction_data = disruptive_prediction_data

    def get_disruptive_interactors(self, protein, mutation):
        disruptive_predicted_interactors = self.disruptive_prediction_data[
            (self.disruptive_prediction_data["UniProt_ID"] == protein) &
            (self.disruptive_prediction_data["Mutation"] == mutation) &
            (self.disruptive_prediction_data["Prediction"] == ClassLabels.DISRUPTING)  # Just to be double sure.
            ]["Interactor_UniProt_ID"].to_list()

        return disruptive_predicted_interactors

    def get_non_disruptive_interactors(self, protein, mutation):
        non_disruptive_predicted_interactors = self.prediction_data[
            (self.prediction_data["UniProt_ID"] == protein) &
            (self.prediction_data["Mutation"] == mutation) &
            (self.prediction_data["Prediction"] == ClassLabels.NONDISRUPTING)  # increasing or no effect (Nondisrupting)
            ]["Interactor_UniProt_ID"].to_list()

        return non_disruptive_predicted_interactors

    def get_all_interactors(self, protein, mutation):
        interactors = self.prediction_data[
            (self.prediction_data["UniProt_ID"] == protein) &
            (self.prediction_data["Mutation"] == mutation)
            ]["Interactor_UniProt_ID"].to_list()

        return interactors

    def find_all_interactions(self):
        log.info("Finding all interactions (disruptive and non-disruptive) for each patient ..")
        patient_to_interactions = {}
        for patient in tqdm(self.patients):
            if self.verbose:
                log.debug(f"\tPATIENT: {patient}")
            interactions = []
            patient_snv_data = self.patient_to_snv_data[patient]

            for index, row in patient_snv_data.iterrows():
                protein = row["SWISSPROT"]
                mutation = row["HGVSp_Short"]
                # for current protein.mutation
                all_interactors = self.get_all_interactors(
                    protein, mutation
                )

                # Add triplets
                for interactor in all_interactors:
                    # log.debug(f"\t\interaction: {protein, mutation, interactor}")

                    interactions.append(
                        (protein, mutation, interactor)
                    )

            patient_to_interactions[patient] = interactions

        self.patient_to_interactions = patient_to_interactions

    def find_disruptive_interactions(self):
        log.info("Finding disruptive interactions for each patient ..")
        patient_to_disruptive_interactions = {}
        for patient in tqdm(self.patients):
            if self.verbose:
                log.debug(f"\tPATIENT: {patient}")
            disruptive_interactions = []
            patient_snv_data = self.patient_to_snv_data[patient]

            for index, row in patient_snv_data.iterrows():
                protein = row["SWISSPROT"]
                mutation = row["HGVSp_Short"]
                # for current protein.mutation
                disruptive_predicted_interactions = self.get_disruptive_interactors(
                    protein, mutation
                )

                # Add disruptive predicted triplets
                for interactor in disruptive_predicted_interactions:
                    # log.debug(f"\t\tDisruptive interaction: {protein, mutation, interactor}")

                    disruptive_interactions.append(
                        (protein, mutation, interactor)
                    )

            patient_to_disruptive_interactions[patient] = disruptive_interactions

        self.patient_to_disruptive_interactions = patient_to_disruptive_interactions

    def print_disruptive_interactions_per_patient(self):
        """
        Prints the disruptive interactions for each patient.
        """
        for patient in self.patients:
            if self.identifier == "uniprot":
                print(f"{patient} -> {self.patient_to_disruptive_interactions[patient]}")

            elif self.identifier == "hugo":
                protein, mutation, interactor = self.patient_to_disruptive_interactions[patient]
                disruptive_interaction_hugo = (
                    f"{protein}:{self.gene_retriever.fetch(protein)}",
                    mutation,
                    f"{interactor}:{self.gene_retriever.fetch(interactor)}",
                )

                print(f"{patient} -> {disruptive_interaction_hugo}")

    def load_uniprot_to_gene_id(self):
        log.info("Loading UniProt ID to Gene ID ..")
        uniprot_to_gene_id = {}

        prediction_data = self.prediction_data.copy()
        self_proteins = list(prediction_data["UniProt_ID"].unique())
        interactor_proteins = list(prediction_data["Interactor_UniProt_ID"].unique())
        proteins = sorted(set(self_proteins + interactor_proteins))

        for protein in tqdm(proteins, desc="Retrieving Gene IDs from UniProt API .. "):
            gene = self.gene_retriever.fetch(protein=protein)
            uniprot_to_gene_id[protein] = gene

        self.uniprot_to_gene_id = uniprot_to_gene_id

        log.info("`uniprot_to_gene_id` loaded. ")

    def is_protein_and_its_mutations_interface_only(
            self,
            protein: str,
            mutations: List[str],
    ):
        protein_interface_only_flag = False
        for mutation in mutations:
            if is_in_elaspic(protein, mutation, self.core_data, self.interface_data):
                if is_core(protein, mutation, self.core_data):
                    return False

                else:
                    protein_interface_only_flag = True

        return protein_interface_only_flag

    def prepare_patient_protein_to_C_I_status_data(
            self,
            patients,
    ):
        log.info(f"Preparing patient_to_CI_status for {len(self.patients)} patients in {self.tcga}.")
        patient_protein_to_C_I_status = dict()

        for patient in tqdm(patients):
            # print(f"Patient: {patient}")
            patient_snv = self.get_patient_snv_data(patient)
            patient_protein_to_mutations = get_patient_protein_to_mutations_dict(patient_snv)

            for protein, mutations in patient_protein_to_mutations.items():

                interface_only = self.is_protein_and_its_mutations_interface_only(
                    protein, mutations
                )

                if interface_only:
                    patient_protein_to_C_I_status[(protein, patient)] = "I"
                    # if protein == "P04637":  # ["Q9Y616", "Q9C0D5", P04637]:
                    #     print(f" - patient adding -- {patient} - {protein}.{mutations}")

                else:
                    patient_protein_to_C_I_status[(protein, patient)] = "C"

        # Prepare patient_protein_mutation_to_C_I_status data
        data = pd.DataFrame(
            patient_protein_to_C_I_status.keys(),
            columns=["PROTEIN", "PATIENT"]
        )
        data["C_I_STATUS"] = patient_protein_to_C_I_status.values()

        self.patient_protein_to_C_I_status_data = data

    def get_current_C_I_status(self, protein, patient):
        c_i_status_data = self.patient_protein_to_C_I_status_data.copy()
        [c_i_status] = c_i_status_data[
            (c_i_status_data["PROTEIN"] == protein) &
            (c_i_status_data["PATIENT"] == patient)
        ]["C_I_STATUS"]

        return c_i_status

    def get_disruptive_prediction_prob(self, protein, mutation, interactor):
        """
        The returned value is the probability of Disruptive class, i.e.
            if ClassLabel.DISRUPTING == 0:
                returns 1 - class_1 (in this case nondisrupting) probability.

            else, that means ClassLabel.DISRUPTING == 1:
                returns class_1 (in this case disrupting) probability

        NOTE:
            Assumed that the Median_Probability is for positive class (class 1) whatever that might be.
        """

        query_data = self.disruptive_prediction_data[
            (self.disruptive_prediction_data["UniProt_ID"] == protein) &
            (self.disruptive_prediction_data["Mutation"] == mutation) &
            (self.disruptive_prediction_data["Interactor_UniProt_ID"] == interactor)
        ]

        if ClassLabels.DISRUPTING == 0:
            [incr_or_no_effect_prediction_probability] = query_data["Median_Probability"]
            # Convert the class nondisruptive probability to class disruptive probability.
            disruptive_prediction_probability = 1 - incr_or_no_effect_prediction_probability
            # Round the number into two decimal places.
            return round(disruptive_prediction_probability, 2)

        elif ClassLabels.DISRUPTING == 1:
            [disruptive_prediction_probability] = query_data["Median_Probability"]
            # Round the number into two decimal places.
            return round(disruptive_prediction_probability, 2)

        else:
            raise ValueError(
                f"Something wrong with ClassLabels. "
                f"Cannot have {ClassLabels.DISRUPTING} of type {type(ClassLabels.DISRUPTING)}"
            )

    def construct_analysis_table(self):
        log.info("Constructing the analysis table ..")
        patient_to_table_entry_set = defaultdict(list)  # Dictionary of list containing unique dictionaries.
        patient_to_seen_protein_mutation_pairs = defaultdict(list)

        for patient in tqdm(self.patients):
            patient_interactions = self.patient_to_interactions[patient]

            for protein, mutation, _ in patient_interactions:
                """
                # unnecessary duplication: but DOES NOT work fine.
                TCGA-A8-A093	P28062	R216W	['P40306']
                TCGA-A8-A093	Q15842	E237K	['Q14654', 'P63252']  
                TCGA-A8-A093	Q15842	E237K	['Q14654', 'P63252']  <- duplicated entries
                """
                # Skipping unnecessary duplication.
                if (protein, mutation) in patient_to_seen_protein_mutation_pairs[patient]:
                    continue

                interactors = self.get_all_interactors(protein, mutation)
                disruptive_interactors = self.get_disruptive_interactors(protein, mutation)
                non_disruptive_interactors = self.get_non_disruptive_interactors(protein, mutation)
                patient_C_I_status = self.get_current_C_I_status(protein=protein, patient=patient)

                print(
                    f"{patient}\t{protein}\t{mutation}\tAll interactors: {interactors}\t{patient_C_I_status}"
                )

                patient_to_seen_protein_mutation_pairs[patient].append((protein, mutation))

                # Adding table info to table_values_set.
                patient_to_table_entry_set[patient].append(
                    {
                        "PATIENT": patient,

                        "PROTEIN_GENE": f"{protein}:{self.gene_retriever.fetch(protein)}",

                        "MUTATION": mutation,

                        # All Interactors
                        "INTERACTORS": ','.join(
                            map(str, map(lambda x: f"{x}:{self.gene_retriever.fetch(x)}", interactors))
                        ),
                        "NUM_INTERACTORS": len(interactors),

                        # Disruptive Interactors
                        "DISRUPTIVE_INTERACTORS": ','.join(
                            map(str,
                                map(
                                    lambda x: f"{x}:"
                                              f"{self.gene_retriever.fetch(x)}:"
                                              f"{self.get_disruptive_prediction_prob(protein, mutation, x)}",
                                    disruptive_interactors)
                                )
                        ),
                        "NUM_DISRUPTIVE_INTERACTORS": len(disruptive_interactors),

                        # Non-Disruptive Interactors
                        "NON_DISRUPTIVE_INTERACTORS": ','.join(
                            map(str, map(lambda x: f"{x}:{self.gene_retriever.fetch(x)}", non_disruptive_interactors))
                        ),
                        "NUM_NON_DISRUPTIVE_INTERACTORS": len(non_disruptive_interactors),

                        # Core+Interface (C) and Interface (I) status
                        "CORE_INTERFACE_VS_INTERFACE_STATUS": patient_C_I_status
                    }
                )

            print('-' * 100)

        analysis_table_entries = []
        for patient, table_entry_set in patient_to_table_entry_set.items():
            for table_entry_dict in table_entry_set:
                # Add table entry dictionary.
                analysis_table_entries.append(table_entry_dict)

        analysis_table = pd.DataFrame(
            analysis_table_entries,
            columns=[
                "PATIENT",
                "PROTEIN_GENE",
                "MUTATION",
                "INTERACTORS",
                "NUM_INTERACTORS",
                "DISRUPTIVE_INTERACTORS",
                "NUM_DISRUPTIVE_INTERACTORS",
                "NON_DISRUPTIVE_INTERACTORS",
                "NUM_NON_DISRUPTIVE_INTERACTORS",
                "CORE_INTERFACE_VS_INTERFACE_STATUS",
            ]
        )

        self.analysis_table = analysis_table
        log.info("Analysis table constructed.")

    def extract(self, folder=None):
        output_file_date = datetime.today().strftime('%Y-%m-%d')
        filename = f"{self.tcga}_patient_interactions_analysis_table_{output_file_date}.xlsx"

        if folder is not None:
            filename = op.join(folder, filename)

        # Ensure the file is not exists before creating to prevent overwriting.
        if op.isfile(filename):
            log.warning(f'File {filename} is already exist.')

        else:
            # Export
            self.analysis_table.to_excel(filename, index=False)
            log.info(f'{filename} is exported.')
