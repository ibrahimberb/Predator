from collections import Counter

from .gene_id_retrieval import GeneIDFetcher
from .loaders import load_snv_datasets

from ..mylogger import get_handler
import logging
from tqdm.auto import tqdm
import pandas as pd

import os.path as op
from datetime import datetime

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class PatientInteractionAnalysis:

    DATA_DIR = "helpers/helpers_analysis/gene_retrieval"

    def __init__(
            self,
            tcga,
            tcga_snv_path,
            patient_interaction_data_path,
            mapping_data_path=None,
    ):
        self.tcga = tcga.upper()
        self.tcga_snv_path = tcga_snv_path
        self.patient_interaction_data_path = patient_interaction_data_path
        self.snv_data = None
        self.patients = None
        self.patient_to_snv_data = None
        self.patient_interaction_data = None

        if mapping_data_path is None:
            self.MAPPING_DATA_NAME = op.join(self.DATA_DIR, "UNIPROT_GENE_MAPPING.csv")
        else:
            self.MAPPING_DATA_NAME = mapping_data_path

        self.gene_id_fetcher = GeneIDFetcher(self.MAPPING_DATA_NAME)

        self.load_snv_data()
        self.load_patient_ids()
        self.load_patient_to_snv_data()
        self.load_patient_interaction_data()

    def load_snv_data(self):
        log.debug("Loading SNV data simplified ..")
        data_materials = {}
        load_snv_datasets(self.tcga, self.tcga_snv_path, data_materials)
        self.snv_data = data_materials[f"{self.tcga}_snv_data_simplified"]

    def load_patient_ids(self):
        log.debug("Loading patient ids ..")
        patients = list(self.snv_data["Tumor_Sample_Barcode"].unique())
        self.patients = patients

    def get_patient_snv_data(self, patient_id):
        patient_snv_data = self.snv_data[
            self.snv_data["Tumor_Sample_Barcode"] == patient_id
            ]
        return patient_snv_data

    def load_patient_to_snv_data(self):
        log.debug("Loading patient to snv_data ..")
        patient_to_snv_data = {}
        for patient in tqdm(self.patients):
            patient_snv_data = self.get_patient_snv_data(patient)
            patient_to_snv_data[patient] = patient_snv_data

        self.patient_to_snv_data = patient_to_snv_data

    def load_patient_interaction_data(self):
        log.debug("patient interaction data patient data ..")
        patient_interaction_data = pd.read_excel(self.patient_interaction_data_path)
        patient_interaction_data = patient_interaction_data[
            patient_interaction_data["CORE_INTERFACE_VS_INTERFACE_STATUS"] == "I"
        ]
        self.patient_interaction_data = patient_interaction_data

    def get_patients_with(self, identifier, identifier_type):

        if identifier_type == "protein":
            patients = self.snv_data[
                self.snv_data["SWISSPROT"] == identifier
                ]["Tumor_Sample_Barcode"].unique()

            return list(patients)

        elif identifier_type == "gene":
            # not that simple. cannot use fetcher because multiple protein id may occur for the same GeneID.
            raise NotImplementedError

        else:
            raise ValueError(f"Invalid argument `{identifier_type}` for `identifier_type`.")

    def calculate_mutual_exclusivity(self, identifier_1, identifier_2, verbose=False, return_num_patients=False):
        """
                             | S1 union S2 |
        Mutual exclusivity = ---------------
                              |S1|  +  |S2|

        """
        # todo: using genes to SNV data may not be a good idea
        #  an alternative is using Protein, and retrieving gene id of that protein, and then querying in SNV data.
        s1 = set(self.get_patients_with(identifier_1, "protein"))
        s2 = set(self.get_patients_with(identifier_2, "protein"))

        if verbose:
            print(f"S1: {s1} ({len(s1)} patients)")
            print(f"S2: {s2} ({len(s2)} patients)")

        mutual_exclusivity_value = len(s1.union(s2)) / (len(s1) + len(s2))

        if return_num_patients:
            return mutual_exclusivity_value, len(s1), len(s2)

        else:
            return mutual_exclusivity_value

    def get_disrupted_interactors_probabilities(self, identifier, identifier_type):
        if identifier_type == "protein":
            identifier_pos = 0
        elif identifier_type == "gene":
            identifier_pos = 1
        else:
            raise ValueError(f"Invalid argument `{identifier_type}` for `identifier_type`.")

        patient_data = self.patient_interaction_data.copy()

        identifier_disruptive_interactors_series = patient_data[
            (patient_data["PROTEIN_GENE"].apply(lambda x: x.split(':')[identifier_pos]) == identifier)
        ]["DISRUPTIVE_INTERACTORS"]

        # Drop nan entries
        identifier_disruptive_interactors_series = identifier_disruptive_interactors_series.dropna()

        # Explode values
        identifier_disruptive_interactors_series = identifier_disruptive_interactors_series.apply(
            lambda x: x.split(",")
        ).explode()

        def _split_prob(x):
            protein_gene, probability = x.rsplit(":", 1)
            probability = float(probability)
            return protein_gene, probability

        identifier_disruptive_interactors_list = list(
            identifier_disruptive_interactors_series.apply(lambda x: _split_prob(x))
        )

        identifier_disruptive_interactors_data = pd.DataFrame(
            identifier_disruptive_interactors_list, columns=["PROTEIN:GENE", "PROBABILITY"]
        )

        identifier_disruptive_interactors_proba_sum_per_interactor = (
            identifier_disruptive_interactors_data.groupby("PROTEIN:GENE")["PROBABILITY"].sum()
        )

        identifier_disruptive_interactors_proba_sum_per_interactor_dict = dict(
            identifier_disruptive_interactors_proba_sum_per_interactor
        )

        # Sort dictionary by value
        identifier_disruptive_interactors_proba_sum_per_interactor_dict_sorted = dict(
            sorted(
                identifier_disruptive_interactors_proba_sum_per_interactor_dict.items(),
                key=lambda item: item[1],
                reverse=True
            )
        )

        return identifier_disruptive_interactors_proba_sum_per_interactor_dict_sorted

    def get_disrupted_interactors(self, identifier, identifier_type, return_counter=False, most_common=False):
        if identifier_type == "protein":
            identifier_pos = 0
        elif identifier_type == "gene":
            identifier_pos = 1
        else:
            raise ValueError(f"Invalid argument `{identifier_type}` for `identifier_type`.")

        patient_data = self.patient_interaction_data.copy()

        identifier_disruptive_interactors_series = patient_data[
            (patient_data["PROTEIN_GENE"].apply(lambda x: x.split(':')[identifier_pos]) == identifier)
        ]["DISRUPTIVE_INTERACTORS"]

        # Drop nan entries
        identifier_disruptive_interactors_series = identifier_disruptive_interactors_series.dropna()

        # Explode values
        identifier_disruptive_interactors_series = identifier_disruptive_interactors_series.apply(
            lambda x: x.split(",")
        ).explode()

        # get rid of probability value.
        identifier_disruptive_interactors = list(
            identifier_disruptive_interactors_series.apply(lambda x: ":".join(x.split(":")[:-1]))
        )

        identifier_disruptive_interactors_counter = Counter(identifier_disruptive_interactors)

        identifier_disruptive_interactors_unique = [
            f"{interactor.split(':')[0]}:{interactor.split(':')[1]}"
            for interactor, _ in identifier_disruptive_interactors_counter.most_common()
        ]

        if return_counter:
            if most_common:
                return identifier_disruptive_interactors_counter.most_common()

            return identifier_disruptive_interactors_counter

        return identifier_disruptive_interactors_unique

    def get_counts_summary_table_protein(self, protein_A):
        interactor_count_pairs = self.get_disrupted_interactors(
            protein_A, identifier_type="protein", return_counter=False
        )

        log.debug(f"interactors: {interactor_count_pairs}")

        entries = []
        for interactor_protein_protein in interactor_count_pairs:
            interactor_protein, _  = interactor_protein_protein.split(':')
            patients_tuple = self.get_patients_with_disruptive_interaction_protein(
                protein_A=protein_A, protein_B=interactor_protein
            )
            patients_protein_A_disrupts_B, patients_protein_B_disrupts_A, patients_intersection = patients_tuple
            entries.append(
                (
                    len(patients_protein_A_disrupts_B),
                    patients_protein_A_disrupts_B,
                    len(patients_protein_B_disrupts_A),
                    patients_protein_B_disrupts_A,
                    len(patients_intersection),
                    patients_intersection
                )
            )

        # Will be concating two dataframes side by side.
        data_1 = pd.DataFrame(
            self.get_disrupted_interactors(
                protein_A, identifier_type="protein", return_counter=True, most_common=True
            ), columns=["PROTEIN_GENE_B", "GENERAL_OCCURRENCE"]
        )
        # Insert TCGA type as the first column value.
        data_1.insert(0, "TCGA", f"{self.tcga}")
        # Insert protein_A as the second column value.
        data_1.insert(1, "PROTEIN_A", f"{protein_A}")
        # Insert gene ID of protein_A as the third column value.
        gene_A_converted = self.gene_id_fetcher.fetch(protein=protein_A)
        data_1.insert(2, "GENE_A", f"{gene_A_converted}")

        data_2 = pd.DataFrame(
            entries, columns=[
                "#_PATIENTS_A_DISR_B",
                "PATIENTS_A_DISR_B",
                "#_PATIENTS_B_DISR_A",
                "PATIENTS_B_DISR_A",
                "#_PATIENTS_INTERSECTION",
                "PATIENTS_INTERSECTION"
            ]
        )

        data_concated = pd.concat([data_1, data_2], axis=1)

        # notify when GENERAL OCCURRENCE differs from number of patients: this happens when
        # a protein is mutated multiple times in the same patient.
        count_difference = data_concated[
            data_concated["GENERAL_OCCURRENCE"] != data_concated["#_PATIENTS_A_DISR_B"]
            ]
        if not count_difference.empty:
            log.warning(
                f"There is a difference between counts."
                f"PROTEIN: {protein_A} INTERACTOR(S): {', '.join(count_difference['PROTEIN_GENE_B'])}"
            )

        # display(data_concated)
        return data_concated

    def get_counts_summary_table(self, gene_A):
        interactor_count_pairs = self.get_disrupted_interactors(
            gene_A, identifier_type="gene", return_counter=False
        )

        log.debug(f"interactors: {interactor_count_pairs}")

        entries = []
        for interactor_protein_gene in interactor_count_pairs:
            _, interactor_gene = interactor_protein_gene.split(':')
            patients_tuple = self.get_patients_with_disruptive_interaction(gene_A=gene_A, gene_B=interactor_gene)
            patients_gene_A_disrupts_B, patients_gene_B_disrupts_A, patients_intersection = patients_tuple
            entries.append(
                (
                    len(patients_gene_A_disrupts_B),
                    patients_gene_A_disrupts_B,
                    len(patients_gene_B_disrupts_A),
                    patients_gene_B_disrupts_A,
                    len(patients_intersection),
                    patients_intersection
                )
            )

        # Will be concating two dataframes side by side.
        data_1 = pd.DataFrame(
            self.get_disrupted_interactors(
                gene_A, identifier_type="gene", return_counter=True, most_common=True
            ), columns=["PROTEIN_GENE_B", "GENERAL_OCCURRENCE"]
        )
        # Insert TCGA type as the first column value.
        data_1.insert(0, "TCGA", f"{self.tcga}")
        # Insert gene_A as the second column value.
        data_1.insert(1, "GENE_A", f"{gene_A}")

        data_2 = pd.DataFrame(
            entries, columns=[
                "#_PATIENTS_A_DISR_B",
                "PATIENTS_A_DISR_B",
                "#_PATIENTS_B_DISR_A",
                "PATIENTS_B_DISR_A",
                "#_PATIENTS_INTERSECTION",
                "PATIENTS_INTERSECTION"
            ]
        )

        data_concated = pd.concat([data_1, data_2], axis=1)

        # notify when GENERAL OCCURRENCE differs from number of patients: this happens when
        # a gene is mutated multiple times in the same patient.
        count_difference = data_concated[
            data_concated["GENERAL_OCCURRENCE"] != data_concated["#_PATIENTS_A_DISR_B"]
        ]
        if not count_difference.empty:
            log.warning(
                f"There is a difference between counts."
                f"GENE: {gene_A} INTERACTOR(S): {', '.join(count_difference['PROTEIN_GENE_B'])}"
            )

        # display(data_concated)
        return data_concated

    def export_counts_summary_table_protein(self, folder_path, protein_A, file_extension="csv"):

        summary_table_protein_A = self.get_counts_summary_table_protein(protein_A)

        gene_A_converted = self.gene_id_fetcher.fetch(protein=protein_A)
        log.debug(f"Exporting Counts Summary Table {self.tcga} {protein_A} {gene_A_converted}..")
        file_name = f"{self.tcga}_{protein_A}_{gene_A_converted}"
        file_name = op.join(folder_path, file_name)
        file_date = datetime.today().strftime('%Y-%m-%d')
        file_name = f'{file_name}_{file_date}.{file_extension}'

        # Ensure the file is not exists before creating to prevent overwriting.
        if op.isfile(file_name):
            log.warning(f"File {file_name} is already exist.\n"
                        "To overwrite existing file, use `overwrite=True`.")
        else:
            # Export
            summary_table_protein_A.to_csv(file_name, index=False)
            log.info(f'{file_name} is exported successfully.')

    def export_counts_summary_table(self, folder_path, gene_A, file_extension="csv"):

        summary_table_gene_A = self.get_counts_summary_table(gene_A)

        log.debug(f"Exporting Counts Summary Table {self.tcga} {gene_A} ..")
        file_name = f"{self.tcga}_{gene_A}"
        file_name = op.join(folder_path, file_name)
        file_date = datetime.today().strftime('%Y-%m-%d')
        file_name = f'{file_name}_{file_date}.{file_extension}'

        # Ensure the file is not exists before creating to prevent overwriting.
        if op.isfile(file_name):
            log.warning(f"File {file_name} is already exist.\n"
                        "To overwrite existing file, use `overwrite=True`.")
        else:
            # Export
            summary_table_gene_A.to_csv(file_name, index=False)
            log.info(f'{file_name} is exported successfully.')

    # This is deprecated ..
    def get_patients_with_disruptive_interaction(self, gene_A, gene_B, verbose=False):
        """
        :param gene_A: Host gene
        :param gene_B: Interactor gene
        :return: patients_gene_A_disrupts_B, patients_gene_B_disrupts_A, patients_intersection
        """
        raise DeprecationWarning
        def _is_gene_exist(x, search_gene) -> bool:
            """
            Checks if any interactions in string x is the search gene.
            :param x:
            :param search_gene:
            :return:
            """
            # Handling `nan` values.
            if isinstance(x, float):
                return False

            protein_gene_prob_list = x.split(',')
            for protein_gene_prob in protein_gene_prob_list:
                _, gene, _ = protein_gene_prob.split(':')
                if search_gene == gene:
                    return True

            return False

        # Patient interaction data with given gene (i.e. filtering for gene_A)
        gene_A_patient_interaction_data = self.patient_interaction_data[
            self.patient_interaction_data["PROTEIN_GENE"].apply(lambda x: x.split(':')[1]) == gene_A
        ].copy()

        # Patient interaction data with given interactor gene. (i.e. filtering for gene_B)
        gene_B_patient_interaction_data = self.patient_interaction_data[
            self.patient_interaction_data["PROTEIN_GENE"].apply(lambda x: x.split(':')[1]) == gene_B
        ].copy()

        # Patients where gene A disrupts gene B
        if gene_A_patient_interaction_data.empty:
            patients_gene_A_disrupts_B = set()
        else:
            patients_gene_A_disrupts_B = set(gene_A_patient_interaction_data[
                gene_A_patient_interaction_data["DISRUPTIVE_INTERACTORS"].apply(lambda x: _is_gene_exist(x, gene_B))
            ]["PATIENT"].unique())

        # Patients where gene B disrupts gene A
        if gene_B_patient_interaction_data.empty:
            patients_gene_B_disrupts_A = set()
        else:
            patients_gene_B_disrupts_A = set(gene_B_patient_interaction_data[
                gene_B_patient_interaction_data["DISRUPTIVE_INTERACTORS"].apply(lambda x: _is_gene_exist(x, gene_A))
            ]["PATIENT"].unique())

        # add intersection.
        patients_intersection = set().intersection(
            patients_gene_A_disrupts_B, patients_gene_B_disrupts_A
        )
        if verbose:
            print(f"GENE_A: {gene_A}")
            print(f"GENE_B: {gene_B}")
            print(f"patients_gene_A_disrupts_B ({len(patients_gene_A_disrupts_B)}): \n {patients_gene_A_disrupts_B}")
            print(f"patients_gene_B_disrupts_A ({len(patients_gene_B_disrupts_A)}): \n {patients_gene_B_disrupts_A}")
            print(f"Intersection ({len(patients_intersection)}): \n{patients_intersection}")
            print('- - - - - - - - - - - - - - - - - - - - - -')

        return (
            patients_gene_A_disrupts_B, patients_gene_B_disrupts_A, patients_intersection
        )

    # todo: this functions can be generalized with `identifier`
    def get_patients_with_disruptive_interaction_protein(self, protein_A, protein_B, verbose=False):
        """
        :param protein_A: Host protein
        :param protein_B: Interactor protein
        :return: patients_protein_A_disrupts_B, patients_protein_B_disrupts_A, patients_intersection
        """

        identifier_pos = 0  # protein

        def _is_protein_exist(x, search_protein) -> bool:
            """
            Checks if any interactions in string x is the search protein.
            :param x:
            :param search_protein:
            :return:
            """
            # Handling `nan` values.
            if isinstance(x, float):
                return False

            protein_gene_prob_list = x.split(',')
            for protein_gene_prob in protein_gene_prob_list:
                protein, _, _ = protein_gene_prob.split(':')
                if search_protein == protein:
                    return True

            return False

        # Patient interaction data with given protein (i.e. filtering for protein_A)
        protein_A_patient_interaction_data = self.patient_interaction_data[
            self.patient_interaction_data["PROTEIN_GENE"].apply(lambda x: x.split(':')[identifier_pos]) == protein_A
            ].copy()

        # Patient interaction data with given interactor protein. (i.e. filtering for protein_B)
        protein_B_patient_interaction_data = self.patient_interaction_data[
            self.patient_interaction_data["PROTEIN_GENE"].apply(lambda x: x.split(':')[identifier_pos]) == protein_B
            ].copy()

        # Patients where protein A disrupts protein B
        if protein_A_patient_interaction_data.empty:
            patients_protein_A_disrupts_B = set()
        else:
            patients_protein_A_disrupts_B = set(
                protein_A_patient_interaction_data[
                    protein_A_patient_interaction_data["DISRUPTIVE_INTERACTORS"].apply(
                        lambda x: _is_protein_exist(x, protein_B)
                    )
                ]["PATIENT"].unique())

        # Patients where protein B disrupts protein A
        if protein_B_patient_interaction_data.empty:
            patients_protein_B_disrupts_A = set()
        else:
            patients_protein_B_disrupts_A = set(
                protein_B_patient_interaction_data[
                    protein_B_patient_interaction_data["DISRUPTIVE_INTERACTORS"].apply(
                        lambda x: _is_protein_exist(x, protein_A)
                    )
                ]["PATIENT"].unique())

        # add intersection.
        patients_intersection = set().intersection(
            patients_protein_A_disrupts_B, patients_protein_B_disrupts_A
        )
        if verbose:
            print(f"PROTEIN_A: {protein_A}")
            print(f"PROTEIN_B: {protein_B}")
            print(f"patients_protein_A_disrupts_B ({len(patients_protein_A_disrupts_B)}): \n {patients_protein_A_disrupts_B}")
            print(f"patients_protein_B_disrupts_A ({len(patients_protein_B_disrupts_A)}): \n {patients_protein_B_disrupts_A}")
            print(f"Intersection ({len(patients_intersection)}): \n{patients_intersection}")
            print('- - - - - - - - - - - - - - - - - - - - - -')

        return (
            patients_protein_A_disrupts_B, patients_protein_B_disrupts_A, patients_intersection
        )

    def get_disruptive_mutual_exclusivity_data(self, protein):
        # uses counts.
        gene = self.gene_id_fetcher.fetch(protein=protein)
        identifier_1_disruptive_interactors_to_values = self.get_disrupted_interactors(
            protein, identifier_type="protein", return_counter=True
        )
        log.info(f"Calculating Mutual Exclusivity over {protein}'s interactors ..")
        log.debug(f"{protein} have {len(identifier_1_disruptive_interactors_to_values)} interactors:\n"
                  f"{identifier_1_disruptive_interactors_to_values}")

        mutual_exclusivity_entries = []
        for interactor, disrupted_interactor_count in identifier_1_disruptive_interactors_to_values.items():
            interactor_protein, interactor_gene = interactor.split(":")
            mut_ex_value, len_s1, len_s2 = self.calculate_mutual_exclusivity(
                protein, interactor_protein, return_num_patients=True
            )

            patients_tuple = self.get_patients_with_disruptive_interaction_protein(
                protein_A=protein, protein_B=interactor_protein
            )
            patients_protein_A_disrupts_B, patients_protein_B_disrupts_A, patients_intersection = patients_tuple

            # log.debug(f"{interactor} \t {mut_ex_value}")
            mutual_exclusivity_entries.append(
                (
                    self.tcga,
                    f"{protein}:{gene}",
                    len_s1,
                    interactor,
                    len_s2,
                    disrupted_interactor_count,
                    len(patients_protein_A_disrupts_B),
                    patients_protein_A_disrupts_B,
                    len(patients_protein_B_disrupts_A),
                    patients_protein_B_disrupts_A,
                    len(patients_intersection),
                    patients_intersection,
                    round(mut_ex_value, 4)),

            )

        mutual_exclusivity_data = pd.DataFrame(
            mutual_exclusivity_entries,
            columns=[
                "TCGA",
                "PROTEIN:GENE",
                "NUM_PATIENTS",
                "INTERACTOR",
                "NUM_PATIENTS_INTERACTOR",
                "DISRUPTIVE_INTERACTOR_COUNT",
                "#_PATIENTS_A_DISR_B",
                "PATIENTS_A_DISR_B",
                "#_PATIENTS_B_DISR_A",
                "PATIENTS_B_DISR_A",
                "#_PATIENTS_INTERSECTION",
                "PATIENTS_INTERSECTION",
                "MUTUAL_EXCLUSIVITY",
            ]
        )

        mutual_exclusivity_data = mutual_exclusivity_data.sort_values("DISRUPTIVE_INTERACTOR_COUNT", ascending=False)

        return mutual_exclusivity_data

    def get_disruptive_mutual_exclusivity_proba_data(self, protein):
        raise DeprecationWarning("not updated as non-proba version..")
        # uses probability sums.
        gene = self.gene_id_fetcher.fetch(protein=protein)
        identifier_1_disruptive_interactors_to_values = self.get_disrupted_interactors_probabilities(
            protein, identifier_type="protein"
        )
        log.info(f"Calculating Mutual Exclusivity over {protein}'s interactors ..")
        log.debug(f"{protein} have {len(identifier_1_disruptive_interactors_to_values)} interactors:\n"
                  f"{identifier_1_disruptive_interactors_to_values}")

        mutual_exclusivity_entries = []
        for interactor, disrupted_interactor_proba_sum in identifier_1_disruptive_interactors_to_values.items():
            interactor_protein, interactor_gene = interactor.split(":")
            mut_ex_value, len_s1, len_s2 = self.calculate_mutual_exclusivity(
                protein, interactor_protein, return_num_patients=True
            )
            # log.debug(f"{interactor} \t {mut_ex_value}")
            mutual_exclusivity_entries.append(
                (
                    f"{protein}:{gene}",
                    len_s1,
                    interactor,
                    len_s2,
                    round(disrupted_interactor_proba_sum, 4),
                    round(mut_ex_value, 4)),

            )

        mutual_exclusivity_data = pd.DataFrame(
            mutual_exclusivity_entries,
            columns=[
                "PROTEIN:GENE",
                "NUM_PATIENTS",
                "INTERACTOR",
                "NUM_PATIENTS_INTERACTOR",
                "DISRUPTIVE_INTERACTOR_PROBA_SUM",
                "MUTUAL_EXCLUSIVITY",
            ]
        )

        mutual_exclusivity_data = mutual_exclusivity_data.sort_values("DISRUPTIVE_INTERACTOR_PROBA_SUM", ascending=False)

        return mutual_exclusivity_data

    def export_disruptive_mutual_exclusivity_data(self, folder_path, protein, file_extension="csv", prob=False):
        # Get Mut Ex data for given identifier.
        if prob:
            mutual_exclusivity_data = self.get_disruptive_mutual_exclusivity_proba_data(protein)
        else:
            mutual_exclusivity_data = self.get_disruptive_mutual_exclusivity_data(protein)

        gene = self.gene_id_fetcher.fetch(protein=protein)

        log.debug(f"Exporting Mutual Exclusivity {self.tcga} {protein} ..")
        file_name = f"{self.tcga}_{protein}_{gene}"
        file_name = op.join(folder_path, file_name)
        file_date = datetime.today().strftime('%Y-%m-%d')

        if prob:
            file_name = f'{file_name}_PROBA'

        file_name = f'{file_name}_{file_date}.{file_extension}'

        # Ensure the file is not exists before creating to prevent overwriting.
        if op.isfile(file_name):
            log.warning(f"File {file_name} is already exist.\n"
                        "To overwrite existing file, use `overwrite=True`.")
        else:
            # Export
            mutual_exclusivity_data.to_csv(file_name, index=False)
            log.info(f'{file_name} is exported successfully.')
