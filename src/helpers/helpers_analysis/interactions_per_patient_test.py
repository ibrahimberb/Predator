import os
import time
from pathlib import Path
from unittest import TestCase
import re
import pandas as pd
import numpy as np
import requests
from tqdm.auto import tqdm

from .interactions_per_patient import InteractionsPerPatient

from ..helpers_analysis.get_protein_to_gene_dict import get_protein_to_gene_dict_via_snv

from ..labels import ClassLabels

from ..mylogger import get_handler
import logging

handler = get_handler(log_type="module")

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class TestInteractionsPerPatient(TestCase):
    # os.chdir("../")

    BRCA_PREDICTION_ID = "ed35a3a3/"
    BRCA_PREDICTIONS_COMMON_PATH = "../data/predictions_datasets/brca_prediction_2022-06-17/" + BRCA_PREDICTION_ID
    PREDICTION_BRCA_REDUCED_PATH = BRCA_PREDICTIONS_COMMON_PATH + "predictions_soft_2022-06-17.csv"
    BRCA_ELASPIC_CORE_PATH = Path("../data/Elaspic_merged_results/BRCA_Core_2021-11-17.txt")
    BRCA_ELASPIC_INTERFACE_PATH = Path("../data/Elaspic_merged_results/BRCA_Interface_2021-11-17.txt")
    PREDICTION_BRCA_REDUCED_PATH = Path(PREDICTION_BRCA_REDUCED_PATH)

    OV_PREDICTION_ID = "865d1897/"
    OV_PREDICTIONS_COMMON_PATH = "../data/predictions_datasets/ov_prediction_2022-06-17/" + OV_PREDICTION_ID
    PREDICTION_OV_REDUCED_PATH = OV_PREDICTIONS_COMMON_PATH + "predictions_soft_2022-06-17.csv"
    OV_ELASPIC_CORE_PATH = Path("../data/Elaspic_merged_results/OV_Core_2021-11-17.txt")
    OV_ELASPIC_INTERFACE_PATH = Path("../data/Elaspic_merged_results/OV_Interface_2021-11-17.txt")
    PREDICTION_OV_REDUCED_PATH = Path(PREDICTION_OV_REDUCED_PATH)

    SNV_COMMON_PATH = "../data/snv_datasets/"
    BRCA_SNV_PATH = Path(os.path.join(SNV_COMMON_PATH, "SNV_BRCA_hg38_2021-09-22.csv"))
    OV_SNV_PATH = Path(os.path.join(SNV_COMMON_PATH, "SNV_OV_hg38_2021-09-22.csv"))

    SNV_PATH = None
    TCGA_PREDICTION_REDUCED_PATH = None
    TCGA_ELASPIC_CORE_PATH = None
    TCGA_ELASPIC_INTERFACE_PATH = None

    TCGA_TEST = "brca"

    def setUp(self) -> None:

        if self.TCGA_TEST == "brca":
            self.SNV_PATH = self.BRCA_SNV_PATH
            self.TCGA_PREDICTION_REDUCED_PATH = self.PREDICTION_BRCA_REDUCED_PATH
            self.TCGA_ELASPIC_CORE_PATH = self.BRCA_ELASPIC_CORE_PATH
            self.TCGA_ELASPIC_INTERFACE_PATH = self.BRCA_ELASPIC_INTERFACE_PATH

        elif self.TCGA_TEST == "ov":
            self.SNV_PATH = self.OV_SNV_PATH
            self.TCGA_PREDICTION_REDUCED_PATH = self.PREDICTION_OV_REDUCED_PATH
            self.TCGA_ELASPIC_CORE_PATH = self.OV_ELASPIC_CORE_PATH
            self.TCGA_ELASPIC_INTERFACE_PATH = self.OV_ELASPIC_INTERFACE_PATH

        else:
            raise

        self.inter_per_pat = InteractionsPerPatient(
            tcga=self.TCGA_TEST,
            prediction_data_path=self.TCGA_PREDICTION_REDUCED_PATH,
            tcga_snv_path=self.SNV_PATH,
            elaspic_core_path=self.TCGA_ELASPIC_CORE_PATH,
            elaspic_interface_path=self.TCGA_ELASPIC_INTERFACE_PATH,
            identifier="uniprot",

        )

    def test_load_snv_data_simplified(self):
        snv_data_simplified = self.inter_per_pat.snv_data_simplified
        columns_list = ["Hugo_Symbol", "SWISSPROT", "HGVSp_Short", "Tumor_Sample_Barcode"]
        self.assertEqual(snv_data_simplified.columns.to_list(), columns_list)

    def test_get_patient_snv_data(self):
        patients = self.inter_per_pat.patients
        for patient in patients:
            patient_snv_data = self.inter_per_pat.get_patient_snv_data(
                patient
            )

            self.assertEqual(patient_snv_data["Tumor_Sample_Barcode"].nunique(), 1)

            self.assertEqual(
                list(patient_snv_data["Tumor_Sample_Barcode"].unique())[0], patient
            )

    def test_patient_ids_are_valid(self):
        patient_ids = self.inter_per_pat.patients
        match = re.compile(r"^TCGA-(\w\w)-(\w\w\w\w)$")
        for patient in patient_ids:
            self.assertIsNotNone(
                match.match(patient), "Test value is not none."
            )

    def test_load_patient_to_snv_data(self):
        patient_to_snv_data = self.inter_per_pat.patient_to_snv_data
        patients = self.inter_per_pat.patients

        self.assertEqual(
            type(patient_to_snv_data), dict
        )

        self.assertEqual(
            sorted(patient_to_snv_data.keys()), sorted(patients)
        )

    def test_load_prediction_data(self):
        prediction_data = self.inter_per_pat.prediction_data

        self.assertEqual(prediction_data.empty, False)

        self.assertEqual(
            set(prediction_data["Prediction"].value_counts().index),
            {ClassLabels.DISRUPTING, ClassLabels.NONDISRUPTING}
        )

    def test_load_disruptive_prediction_data(self):
        disruptive_prediction_data = self.inter_per_pat.disruptive_prediction_data

        self.assertEqual(
            set(disruptive_prediction_data["Prediction"].value_counts().index),
            {ClassLabels.DISRUPTING}
        )

    def test_get_disruptive_predicted_interactions(self):

        if self.inter_per_pat.tcga == "BRCA":

            if self.inter_per_pat.identifier == "uniprot":
                self.assertNotEqual(
                    self.inter_per_pat.get_disruptive_interactors(
                        "Q01196", "G95R"
                    ), []
                )
                self.assertNotEqual(
                    self.inter_per_pat.get_disruptive_interactors(
                        "P00747", "D665H"
                    ), []
                )

            elif self.inter_per_pat.identifier == "hugo":
                pass

        elif self.inter_per_pat.tcga == "OV":

            if self.inter_per_pat.identifier == "uniprot":
                self.assertNotEqual(
                    self.inter_per_pat.get_disruptive_interactors(
                        "O75175", "Q684H"
                    ), []
                )
                self.assertNotEqual(
                    self.inter_per_pat.get_disruptive_interactors(
                        "Q14814", "K31T"
                    ), []
                )

        else:
            raise ValueError("TCGA falls outside test cases.")

        # a dummy fail case.
        self.assertEqual(
            self.inter_per_pat.get_disruptive_interactors(
                "P123123", "MUT123"
            ), []
        )

    # FIXME -- function name is not correct
    def test_patient_to_disruptive_interactions(self):
        patient_to_disruptive_interactions = self.inter_per_pat.patient_to_disruptive_interactions

        if self.inter_per_pat.tcga == "BRCA":
            self.assertEqual(
                len(patient_to_disruptive_interactions), 985
            )

        elif self.inter_per_pat.tcga == "OV":
            self.assertEqual(
                len(patient_to_disruptive_interactions), 436
            )

        else:
            raise ValueError("TCGA falls outside test cases.")

    def test_find_disruptive_interactions(self):

        patient_to_disruptive_interactions = self.inter_per_pat.patient_to_disruptive_interactions
        prediction_data = self.inter_per_pat.prediction_data
        disruptive_prediction_data = self.inter_per_pat.disruptive_prediction_data

        for patient, disruptive_triplets in patient_to_disruptive_interactions.items():
            patient_snv_data = self.inter_per_pat.get_patient_snv_data(patient)
            for disruptive_triplet in disruptive_triplets:
                protein, mutation, interactor = disruptive_triplet

                patient_search_data = patient_snv_data[
                    (patient_snv_data["SWISSPROT"] == protein) &
                    (patient_snv_data["HGVSp_Short"] == mutation)
                ]

                self.assertFalse(patient_search_data.empty, msg="data is not empty")

                prediction_search_data = prediction_data[
                    (prediction_data["UniProt_ID"] == protein) &
                    (prediction_data["Mutation"] == mutation) &
                    (prediction_data["Interactor_UniProt_ID"] == interactor)
                ]

                # Assert there is only one row.
                self.assertEqual(len(prediction_search_data), 1)

                # Ensure there is only one value in `Prediction` column and get it.
                self.assertEqual(prediction_search_data["Prediction"].nunique(), 1)
                [val] = prediction_search_data["Prediction"].values

                # Assert the prediction is "disruptive"
                self.assertEqual(val, ClassLabels.DISRUPTING)

                disruptive_prediction_search_data = disruptive_prediction_data[
                    (disruptive_prediction_data["UniProt_ID"] == protein) &
                    (disruptive_prediction_data["Mutation"] == mutation) &
                    (disruptive_prediction_data["Interactor_UniProt_ID"] == interactor)
                ]

                # Assert there is only one row.
                self.assertEqual(len(disruptive_prediction_search_data), 1)

                # Ensure there is only one value in `Prediction` column and get it.
                self.assertEqual(disruptive_prediction_search_data["Prediction"].nunique(), 1)
                [val] = disruptive_prediction_search_data["Prediction"].values

                # Assert the prediction is "disruptive"
                self.assertEqual(val, ClassLabels.DISRUPTING)

                # Assert the search datasets are identical.
                self.assertTrue(
                    prediction_search_data.equals(disruptive_prediction_search_data),
                    "datasets are not identical."
                )

    def validate_disruptive_interactions_single_patient(self, patient):
        snv_data = self.inter_per_pat.snv_data_simplified
        disruptive_prediction_data = self.inter_per_pat.disruptive_prediction_data
        patient_snv_data = snv_data[snv_data["Tumor_Sample_Barcode"] == patient]

        assert patient_snv_data.equals(
            self.inter_per_pat.patient_to_snv_data[patient]
        )

        patient_disruptive_interactions = (
            self.inter_per_pat.patient_to_disruptive_interactions[patient]
        )

        for disruptive_interactions in patient_disruptive_interactions:
            protein, mutation, interactor = disruptive_interactions

            self.assertFalse(
                patient_snv_data[
                    (patient_snv_data["SWISSPROT"] == protein) &
                    (patient_snv_data["HGVSp_Short"] == mutation)
                    ].empty,
                msg="patient_snv_data is empty!"
            )

            self.assertFalse(
                disruptive_prediction_data[
                    (disruptive_prediction_data["UniProt_ID"] == protein) &
                    (disruptive_prediction_data["Mutation"] == mutation) &
                    (disruptive_prediction_data["Interactor_UniProt_ID"] == interactor)
                    ].empty,
                msg="patient_snv_data is empty!"
            )

    def test_validate_disruptive_interactions_per_patient(self):
        for patient in self.inter_per_pat.patients:
            self.validate_disruptive_interactions_single_patient(patient)

    def validate_all_interactions_single_patient(self, patient):
        snv_data = self.inter_per_pat.snv_data_simplified
        prediction_data = self.inter_per_pat.prediction_data
        patient_snv_data = snv_data[snv_data["Tumor_Sample_Barcode"] == patient]

        assert patient_snv_data.equals(
            self.inter_per_pat.patient_to_snv_data[patient]
        )

        patient_all_interactions = (
            self.inter_per_pat.patient_to_interactions[patient]
        )

        for interactions in patient_all_interactions:
            protein, mutation, interactor = interactions

            self.assertFalse(
                patient_snv_data[
                    (patient_snv_data["SWISSPROT"] == protein) &
                    (patient_snv_data["HGVSp_Short"] == mutation)
                    ].empty,
                msg="patient_snv_data is empty!"
            )

            self.assertFalse(
                prediction_data[
                    (prediction_data["UniProt_ID"] == protein) &
                    (prediction_data["Mutation"] == mutation) &
                    (prediction_data["Interactor_UniProt_ID"] == interactor)
                    ].empty,
                msg="patient_snv_data is empty!"
            )

    def test_validate_all_interactions_per_patient(self):
        for patient in self.inter_per_pat.patients:
            self.validate_all_interactions_single_patient(patient)

    def test_uniprot_gene_conversion_with_SNV(self):
        """
        The SNV data also contains a column for Hugo Gene ID.
        Here, we test whether we could retrieve correct gene id from Uniprot API
        by comparing it with the one in SNV data.
        """

        log.debug(f"Loading SNV data: {self.SNV_PATH}")
        snv_data = pd.read_csv(self.SNV_PATH, low_memory=False)
        log.debug("SNV data loaded.")

        prediction_data = self.inter_per_pat.prediction_data
        self_proteins = list(prediction_data["UniProt_ID"].unique())
        interactor_proteins = list(prediction_data["Interactor_UniProt_ID"].unique())
        proteins = sorted(set(self_proteins + interactor_proteins))

        log.debug(f"len self_proteins: {len(self_proteins)}")
        log.debug(f"len interactor_proteins: {len(interactor_proteins)}")
        log.debug(f"len proteins: {len(proteins)}")

        protein_to_gene_dict = get_protein_to_gene_dict_via_snv(proteins, snv_data)
        tcga_exception_proteins = []

        for protein in tqdm(proteins):
            log.info(f"\nCURRENT PROTEIN: {protein}")
            if protein_to_gene_dict[protein] == "NA":
                log.error(f"PROTEIN: {protein} IS `NA` IN SNV DATA.")
                converted_gene = self.inter_per_pat.gene_retriever.fetch(protein)
                converted_legacy_gene = self.get_gene_id_from_uniprot_regex_v1(protein)
                log.debug(f"converted_gene:         {converted_gene}")
                log.debug(f"converted_legacy_gene:  {converted_legacy_gene}")

                if protein in ["Q5U077"]:
                    pass

                else:
                    self.assertTrue(
                        (converted_gene == converted_legacy_gene) or
                        (converted_gene is converted_legacy_gene)
                    )

                log.warning(f"Retrieved from UNIPROT API: {converted_gene}")
                continue

            # Known exceptions
            # TODO: extract these into a text file, because we might have more than 100 of these ..
            #  and uses `tcga` too

            # exception_proteins = tcga_exception_proteins

            # if protein in exception_proteins:
            #     log.warning(f"EXCEPTION PROTEIN: {protein}")
            #     converted_gene = self.disruptive_interactions_per_patient.get_gene_id_from_uniprot(protein)
            #     converted_legacy_gene = self.get_gene_id_from_uniprot_regex_v1(protein)
            #     log.debug(f"converted_gene:         {converted_gene}")
            #     log.debug(f"converted_legacy_gene:  {converted_legacy_gene}")
            #     continue

            actual_gene = protein_to_gene_dict[protein]
            converted_gene = self.inter_per_pat.gene_retriever.fetch(protein)
            converted_legacy_gene = self.get_gene_id_from_uniprot_regex_v1(protein)

            log.debug(f"actual_gene:            {actual_gene}")
            log.debug(f"converted_gene:         {converted_gene}")
            log.debug(f"converted_legacy_gene:  {converted_legacy_gene}")

            if converted_gene == "N/A" or converted_legacy_gene == "N/A":
                log.error("CANNOT FIND CONVERTED GENE ID IN DATABASE ..")
                continue

            try:
                self.assertTrue(converted_gene == converted_legacy_gene == actual_gene)
                self.assertEqual(converted_gene, converted_legacy_gene)
                self.assertEqual(actual_gene, converted_gene)
            except AssertionError:
                log.error(f"EXCEPTIONAL PROTEIN: {protein}. ADDING TO LIST ..")
                tcga_exception_proteins.append((protein, actual_gene, converted_gene))

        with open(f"{self.inter_per_pat.tcga}_exceptions.txt", "w") as file:
            log.info("WRITING TO FILE ...")
            for exception_protein, actual_gene, converted_gene in tcga_exception_proteins:
                file.write(f"{exception_protein}\t{actual_gene}\t{converted_gene}\n")

            log.info("WRITING TO FILE COMPLETED.")

        log.info("TEST `test_uniprot_gene_conversion` COMPLETE ...")

    @staticmethod
    def get_gene_id_from_uniprot_regex_v1(uniprot_id):

        def get_gene_from_fasta_regex_v1(fasta_text):
            info_line = fasta_text.split('\n')[0]
            line_splitted = info_line.split()

            pattern = re.compile(r"^GN=(.+)")
            try:
                [gene] = filter(pattern.match, line_splitted)
            except ValueError:
                log.warning(f"NO GENE ID FOUND FOR {uniprot_id}.")
                return np.nan

            gene = gene.replace("GN=", "")

            return gene

        log.debug("Retrieving sequence {} ...".format(uniprot_id))
        address = "http://www.uniprot.org/uniprot/{}.fasta".format(uniprot_id)

        n_attempt = 3
        attempt = 0
        while attempt < n_attempt:
            r = requests.get(address)
            if r.status_code == 200:
                gene = get_gene_from_fasta_regex_v1(r.text)
                return gene

            attempt += 1
            log.warning(f"attempt: {attempt}")
            log.warning(f"status_code: {r.status_code}")
            time.sleep(1)

        log.error(f"COULD NOT RETRIEVE GENE: {attempt}")
        return np.nan
