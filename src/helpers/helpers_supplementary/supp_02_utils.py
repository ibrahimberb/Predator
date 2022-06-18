from typing import List

from pathlib import Path
import pandas as pd

from tqdm.notebook import tqdm

from ..helpers_analysis.loaders import load_snv_datasets, load_elaspic_datasets

from ..mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


def merge_tcga_supp_02_datasets(*datasets):
    merged_supp_02_data = pd.concat(datasets)
    merged_supp_02_data.index.name = "TCGA"
    return merged_supp_02_data


class Supp02Helper:
    def __init__(
            self,
            tcga: str,
            tcga_snv_path: Path,
            elaspic_core_path: Path,
            elaspic_interface_path: Path,
    ):
        self.tcga = tcga.upper()
        self.tcga_snv_path = tcga_snv_path
        self.elaspic_core_path = elaspic_core_path
        self.elaspic_interface_path = elaspic_interface_path

        self.snv_data = None
        self.snv_data_simplified = None
        self.patients = None
        self.patient_to_snv_data = None
        self.elaspic_data_materials = None
        self.load_materials()

        self.analysis_info = {}
        self.analysis_table = None
        self.prepare_info()

    def load_materials(self):
        log.info("Loading materials ..")
        self.load_snv_data()
        self.load_patient_ids()
        self.load_patient_to_snv_data()
        self.load_elaspic_core_and_interface_datasets()
        log.info("Materials loaded.")
        log.info(f"Number of {self.tcga} patients: {len(self.patients)}.")

    def load_snv_data(self):
        log.debug("Loading SNV data simplified ..")
        data_materials = {}
        load_snv_datasets(self.tcga, self.tcga_snv_path, data_materials)
        self.snv_data = data_materials[f"{self.tcga}_snv_data"]
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

    def load_elaspic_core_and_interface_datasets(self):
        elaspic_data_materials = {}
        load_elaspic_datasets(
            tcga=self.tcga,
            elaspic_core_path=self.elaspic_core_path,
            elaspic_interface_path=self.elaspic_interface_path,
            data_materials=elaspic_data_materials,
        )
        self.elaspic_data_materials = elaspic_data_materials

    def prepare_info(self):
        """
        Get number of entries and number of mutations.

        Notes
        -----
            Number of mutations is found based on unique PROTEIN.MUTATIONS in corresponding data.
            When implementing, I use data[~data.duplicated(subset=['PROTEIN', 'MUTATION'])] pattern.
            Notice the '~' symbol, so it is taking non duplicated entries of columns PROTEIN and MUTATION.

        Returns
        -------
            None. Prepares the analysis table.
        """

        num_patients = len(self.patients)

        num_original_snv_entries = len(self.snv_data)
        num_original_snv_mutations = len(
            self.snv_data[~self.snv_data.duplicated(subset=["SWISSPROT", "HGVSp_Short"])]
        )

        num_snv_data_simplified_entries = len(self.snv_data_simplified)
        num_snv_data_simplified_mutations = len(
            self.snv_data_simplified[~self.snv_data_simplified.duplicated(subset=["SWISSPROT", "HGVSp_Short"])]
        )

        # ELASPIC CORE ENTRIES AND MUTATIONS
        temp_df = self.elaspic_data_materials[f"{self.tcga}_elaspic_core_data"]
        num_elaspic_core_entries = len(temp_df)
        num_elaspic_core_mutations = len(
            temp_df[~temp_df.duplicated(subset=["UniProt_ID", "Mutation"])]
        )

        # ELASPIC CORE PROCESSED ENTRIES AND MUTATIONS
        temp_df = self.elaspic_data_materials[f"{self.tcga}_elaspic_core_data_simplified"]
        num_elaspic_core_processed_entries = len(temp_df)
        # most likely as the same as `num_elaspic_core_mutations`.
        num_elaspic_core_processed_mutations = len(
            temp_df[~temp_df.duplicated(subset=["UniProt_ID", "Mutation"])]
        )

        # ELASPIC INTERFACE ENTRIES AND MUTATIONS
        temp_df = self.elaspic_data_materials[f"{self.tcga}_elaspic_interface_data"]
        num_elaspic_interface_entries = len(temp_df)
        num_elaspic_interface_mutations = len(
            temp_df[~temp_df.duplicated(subset=["UniProt_ID", "Mutation"])]
        )

        # ELASPIC INTERFACE PROCESSED TRIPLETS (also entries) AND MUTATIONS
        temp_df = self.elaspic_data_materials[f"{self.tcga}_elaspic_interface_processed_data"]
        num_elaspic_interface_processed_triplets = len(temp_df)
        # most likely as the same as `num_elaspic_interface_mutations`.
        num_elaspic_interface_processed_mutations = len(
            temp_df[~temp_df.duplicated(subset=["UniProt_ID", "Mutation"])]
        )

        # ELASPIC CORE PROCESSED AND INTERFACE PROCESSED ENTRIES AND MUTATIONS
        temp_df = self.elaspic_data_materials[f"{self.tcga}_elaspic_core_and_interface_data"]
        num_elaspic_core_processed_and_interface_processed_entries = len(temp_df)
        num_elaspic_core_processed_and_interface_processed_mutations = len(
            temp_df[~temp_df.duplicated(subset=["UniProt_ID", "Mutation"])]
        )

        analysis_info = {
            "num_patients": num_patients,
            "num_original_snv_entries": num_original_snv_entries,
            "num_original_snv_mutations": num_original_snv_mutations,

            "num_snv_data_simplified_entries": num_snv_data_simplified_entries,
            "num_snv_data_simplified_mutations": num_snv_data_simplified_mutations,

            "num_elaspic_core_entries": num_elaspic_core_entries,
            "num_elaspic_core_mutations": num_elaspic_core_mutations,

            "num_elaspic_core_processed_entries": num_elaspic_core_processed_entries,
            "num_elaspic_core_processed_mutations": num_elaspic_core_processed_mutations,

            "num_elaspic_interface_entries": num_elaspic_interface_entries,
            "num_elaspic_interface_mutations": num_elaspic_interface_mutations,

            "num_elaspic_interface_processed_triplets": num_elaspic_interface_processed_triplets,
            "num_elaspic_interface_processed_mutations": num_elaspic_interface_processed_mutations,

            "num_elaspic_core_processed_and_interface_processed_entries":
                num_elaspic_core_processed_and_interface_processed_entries,
            "num_elaspic_core_processed_and_interface_processed_mutations":
                num_elaspic_core_processed_and_interface_processed_mutations,
        }

        self.analysis_info = analysis_info

        self.analysis_table = pd.DataFrame(self.analysis_info, [f"{self.tcga}"])

        log.info(f"Analysis table prepared for {self.tcga}.")


class Supp02:
    def __init__(self, tcga_list: List[str]):
        pass
