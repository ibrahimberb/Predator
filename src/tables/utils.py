import os

from pandas import DataFrame
import re
import pandas as pd
import os.path as op

from IPython.display import display


def extract_filtered_proteins(data, path):
    [tcga] = data["TCGA"].unique()
    proteins = data["PROTEIN"]
    genes = data["GENE"]
    file_name = f"{tcga}_filtered_proteins.txt"
    file_path = op.join(path, file_name)
    with open(file_path, "w") as file:
        for protein, gene in zip(proteins, genes):
            file.write(f"{protein}_{gene}\n")

    print(f"{len(proteins)} proteins extracted into {file_name}")


class TableIsEmptyError(Exception):
    pass


class MostFrequentlyDisruptedPartnersAdder:
    def __init__(
            self,
            tcga: str,
            preliminary_data: DataFrame,
            interactions_summary_folder_path: str
    ):
        self.tcga = tcga.upper()
        self.preliminary_data = preliminary_data
        self.interactions_summary_folder_path = op.join(
            interactions_summary_folder_path, self.tcga
        )
        self.proteins = self.preliminary_data["PROTEIN"].tolist()

    def get_most_frequently_disrupted_interaction(self, protein):
        try:
            protein_summary_table = self._get_protein_interaction_summary_table(protein)
        except TableIsEmptyError:
            return "NA"

        highest_frequency = protein_summary_table["#_PATIENTS_A_DISR_B"].max()
        protein_summary_table_highest = protein_summary_table[
            protein_summary_table["#_PATIENTS_A_DISR_B"] == highest_frequency
            ].copy()

        # display(protein_summary_table_highest)
        most_frequently_disrupted_partners = list(protein_summary_table_highest["PROTEIN_GENE_B"])
        most_frequently_disrupted_partners = [partner.split(":")[1] for partner in most_frequently_disrupted_partners]
        most_frequently_disrupted_partners = sorted(set(most_frequently_disrupted_partners))

        # If patient count is 1 or 2, we exclude them.
        if highest_frequency <= 2:
            return "-"

        else:
            return f'{", ".join(most_frequently_disrupted_partners)} ({highest_frequency})'

    def get_most_frequently_disrupted_partners_added_data(self):
        most_frequently_disrupted_partners = []
        preliminary_data_updated = self.preliminary_data.copy()
        for protein in self.proteins:
            most_frequently_disrupted_partners.append(
                self.get_most_frequently_disrupted_interaction(protein)
            )

        preliminary_data_updated.insert(
            loc=3,
            column="MOST_FREQUENTLY_DISRUPTED_PARTNER",
            value=most_frequently_disrupted_partners,
        )

        return preliminary_data_updated

    def _get_protein_interaction_summary_table(self, protein):
        files = os.listdir(self.interactions_summary_folder_path)
        files_found = []
        for file in files:
            if re.search(r"{}_{}_.*".format(self.tcga, protein), file):
                files_found.append(file)

        try:
            [file_found] = files_found
            file_path = op.join(self.interactions_summary_folder_path, file_found)
            protein_interaction_summary_table = pd.read_csv(file_path)
            if protein_interaction_summary_table.empty:
                return TableIsEmptyError

            return protein_interaction_summary_table

        except ValueError:
            print(f"files: {files}")
            print(f"protein: {protein}")
            print(f"interactions_summary_folder_path: {self.interactions_summary_folder_path}")
            raise


