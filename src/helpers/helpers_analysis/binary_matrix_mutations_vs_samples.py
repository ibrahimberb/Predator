from datetime import datetime

from .loaders import (
    load_snv_datasets,
)

from pathlib import Path
import os.path as op

import pandas as pd

from ..mylogger import get_handler
import logging

from tqdm.auto import tqdm

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class BinaryMatrix:
    """
    The class for constructing binary matrix that includes mutations versus samples.
    """

    def __init__(
            self,
            tcga: str,
            tcga_snv_path: Path,
            filtering_pairs_filepath=None
    ):
        self.tcga = tcga.upper()
        self.tcga_snv_path = tcga_snv_path
        self.filtering_pairs_filepath = filtering_pairs_filepath

        self.snv_data_simplified = None
        self.patients = None
        self.proteins = None
        self.patient_to_snv_data = None

        self.binary_matrix = None

        self.load_snv_data_simplified()
        self.load_patient_ids()
        self.load_proteins()
        self.load_patient_to_snv_data()

        if self.filtering_pairs_filepath is None:
            self.filtering_proteins = None
        else:
            self.filtering_proteins = self.get_existing_proteins()

        self.construct_binary_matrix()

        if self.filtering_pairs_filepath is None:
            self.test_binary_matrix()

    def get_existing_proteins(self):
        pairs_data = pd.read_csv(
            self.filtering_pairs_filepath,
            sep="_",
            header=None,
            names=["TCGA", "PROTEIN_GENE", "INTERACTOR_PROTEIN_GENE"],
        )

        pairs_data["PROTEIN"] = pairs_data["PROTEIN_GENE"].apply(lambda x: x.split(":")[0])
        pairs_data["INTERACTOR_PROTEIN"] = pairs_data["INTERACTOR_PROTEIN_GENE"].apply(lambda x: x.split(":")[0])

        proteins_exist = sorted(
            set(pairs_data["PROTEIN"]).union(set(pairs_data["INTERACTOR_PROTEIN"]))
        )

        return proteins_exist

    def load_snv_data_simplified(self):
        log.debug("Loading SNV data simplified ..")
        data_materials = {}
        load_snv_datasets(self.tcga, self.tcga_snv_path, data_materials)
        self.snv_data_simplified = data_materials[f"{self.tcga}_snv_data_simplified"]

    def load_patient_ids(self):
        log.debug("Loading patient ids ..")
        patients = sorted(self.snv_data_simplified["Tumor_Sample_Barcode"].unique())
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

    def load_proteins(self):
        log.debug("Loading proteins ..")
        proteins = sorted(self.snv_data_simplified["SWISSPROT"].unique())
        self.proteins = proteins

    def _get_protein_vector_patient(self, patient):
        patient_snv = self.patient_to_snv_data[patient]
        patient_proteins = sorted(patient_snv["SWISSPROT"].unique())
        protein_occurrence_vector = [
            int(protein in patient_proteins) for protein in self.proteins
        ]
        return protein_occurrence_vector

    def construct_binary_matrix(self):
        log.debug("Constructing binary matrix ..")
        binary_matrix_data = pd.DataFrame(
            0, columns=self.patients, index=self.proteins
        )

        for patient in tqdm(self.patients):
            binary_matrix_data[patient] = self._get_protein_vector_patient(patient)

        binary_matrix_data.index.name = "PROTEIN"

        if self.filtering_pairs_filepath is not None:
            overlapping_proteins = sorted(
                set(self.filtering_proteins).intersection(
                    set(binary_matrix_data.index)
                )
            )
            binary_matrix_data = binary_matrix_data.loc[overlapping_proteins, :]

        self.binary_matrix = binary_matrix_data
        log.info("Constructing binary matrix completed.")

    def test_binary_matrix(self):
        for patient in tqdm(self.patients, desc="Running the test .."):
            proteins_in_matrix = self.binary_matrix[
                self.binary_matrix[patient] == 1
            ].index

            proteins_in_snv = self.patient_to_snv_data[patient]["SWISSPROT"].unique()

            assert sorted(proteins_in_matrix) == sorted(proteins_in_snv)

        log.info("All tests pass.")

    def extract(self, folder=None, compress=False):
        output_file_date = datetime.today().strftime('%Y-%m-%d')
        filtered = "_filtered" if self.filtering_pairs_filepath is not None else ""
        filename = f"{self.tcga}_binary_matrix_{output_file_date}{filtered}.csv"

        if folder is not None:
            filename = op.join(folder, filename)

        # Ensure the file is not exists before creating to prevent overwriting.
        if op.isfile(filename):
            log.warning(f'File {filename} is already exist.')

        else:
            # Export
            if compress:
                filename = f"{filename}.gz"
                self.binary_matrix.to_csv(filename, compression="gzip")
                log.info(f'{filename} is exported (compressed).')
            else:
                self.binary_matrix.to_csv(filename)
                log.info(f'{filename} is exported.')
