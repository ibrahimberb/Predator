from pathlib import Path
from typing import List

import pandas as pd

from helpers.helpers_analysis.loaders import (
    load_snv_datasets,
    load_prediction_dataset,
    load_elaspic_datasets,
    load_reference_dataset,
    load_patient_interaction_dataset,
    ReferenceDataset,
)

from helpers.helpers_analysis.get_elaspic_proteins import get_elaspic_proteins
from helpers.helpers_analysis.get_protein_to_gene_dict import get_protein_to_gene_dict
from helpers.helpers_analysis.get_protein_to_num_elaspic_entries_dict import \
    get_protein_to_num_elaspic_interface_entires_dict
from helpers.helpers_analysis.get_protein_to_num_unique_interactors import get_protein_to_num_unique_interactors
from helpers.helpers_analysis.column_adders import (
    add_baseline,
    add_elaspic_coverage,
    add_genes,
    add_num_disruptive_entries,
    add_num_elaspic_interface_entries,
    add_num_incr_noeff_entries,
    add_num_unique_interactors,
    add_our_method,
    add_patient_core_count,
    add_patient_interface_count,
    add_cancermine_status,
    add_num_interface_patients_disruptive_interactor,
    add_cgc_status
)

from helpers.helpers_analysis.counts_baseline_vs_our_method import counts_baseline_vs_our_method
from helpers.helpers_analysis.plot_roc_curve import roc_curve_analysis
from helpers.helpers_analysis.common import save_auc_scores


from helpers.mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)

TCGA_CODE = str
REFERENCE_DATA_NAME = str
RANDOM_ID = str
CohortSpecificReferenceDataPath, ReferenceDataPath = Path, Path


class PredatorAnalysis:
    def __init__(
            self,
            tcga: TCGA_CODE,
            snv_path: Path,
            prediction_data_path: Path,
            prediction_id: RANDOM_ID,
            elaspic_core_path: Path,
            elaspic_interface_path: Path,
            patient_interaction_data_path: Path,
            reference_data_name: REFERENCE_DATA_NAME,
            reference_data_spec_cohort_path: Path,
            reference_data_path: Path,
    ):

        log.info(" - - Predator Analysis - - ")
        log.info(f"TCGA: {tcga}")
        log.info(f"PREDICTION ID: {prediction_id}")

        self.tcga = tcga
        self.snv_path = snv_path
        self.prediction_data_path = prediction_data_path
        self.prediction_id = prediction_id
        self.elaspic_core_path = elaspic_core_path
        self.elaspic_interface_path = elaspic_interface_path
        self.reference_data_name = reference_data_name
        self.reference_data_spec_cohort_path = reference_data_spec_cohort_path
        self.reference_data_path = reference_data_path
        self.patient_interaction_data_path = patient_interaction_data_path
        self.data_materials = {}
        self.auc_scores = {}

        # setattr(self, self.tcga, {})

        self.load_datasets()

        log.info("Initialization completed.\n")

    def load_datasets(self):
        load_snv_datasets(self.tcga, self.snv_path, self.data_materials)
        load_prediction_dataset(self.tcga, self.prediction_data_path, self.data_materials)
        load_elaspic_datasets(
            self.tcga, self.elaspic_core_path, self.elaspic_interface_path, self.data_materials
        )
        load_reference_dataset(
            self.tcga,
            self.reference_data_name,
            self.reference_data_spec_cohort_path,
            self.reference_data_path,
            self.data_materials
        )
        load_patient_interaction_dataset(self.patient_interaction_data_path,  self.data_materials)

    def prepare_analysis(self):
        # Proteins
        elaspic_proteins = get_elaspic_proteins(
            self.data_materials[f"{self.tcga}_elaspic_core_data"],
            self.data_materials[f"{self.tcga}_elaspic_interface_data"]
        )
        self.data_materials[f"{self.tcga}_elaspic_proteins"] = elaspic_proteins
        log.debug(f"{self.tcga}_elaspic_proteins loaded.")
        num_proteins = len(self.data_materials[f"{self.tcga}_elaspic_proteins"])
        log.debug(f"Number of proteins in ELASPIC {self.tcga}: {num_proteins}")

        # Genes
        protein_to_gene_dict = get_protein_to_gene_dict(
            self.data_materials[f"{self.tcga}_elaspic_proteins"]
        )
        self.data_materials[f"{self.tcga}_protein_to_gene_dict"] = protein_to_gene_dict
        log.debug(f"{self.tcga}_protein_to_gene_dict loaded.")

        # ELASPIC Number of Interface Entries
        protein_to_num_elaspic_interface_entries = get_protein_to_num_elaspic_interface_entires_dict(
            self.data_materials[f"{self.tcga}_elaspic_proteins"],
            self.data_materials[f"{self.tcga}_elaspic_interface_processed_data"]
        )
        self.data_materials[f"{self.tcga}_protein_to_num_elaspic_interface_entries"] = protein_to_num_elaspic_interface_entries
        log.debug(f"{self.tcga}_protein_to_num_elaspic_interface_entries loaded.")

        # ELASPIC Number of Unique Interactors
        protein_to_num_unique_interactors = get_protein_to_num_unique_interactors(
            self.data_materials[f"{self.tcga}_elaspic_proteins"],
            self.data_materials[f"{self.tcga}_elaspic_interface_processed_data"]
        )
        self.data_materials[f"{self.tcga}_protein_to_num_unique_interactors"] = protein_to_num_unique_interactors
        log.debug(f"{self.tcga}_protein_to_num_unique_interactors loaded.")

        # Patients
        patients = list(
            self.data_materials[f"{self.tcga}_snv_data_simplified"]['Tumor_Sample_Barcode'].unique()
        )
        log.debug(f'Number of patients in {self.tcga}: {len(patients)}.')
        self.data_materials[f"{self.tcga}_patients"] = patients

    def construct_analysis_table(self):
        # 1. Adding `PROTEIN` Column
        log.debug(f"Adding `PROTEIN` column ..")
        preliminary_data = pd.DataFrame(
            self.data_materials[f"{self.tcga}_elaspic_proteins"], columns=['PROTEIN']
        )

        # 2. Adding `GENE` Column
        log.debug(f"Adding `GENE` column ..")
        add_genes(preliminary_data, self.data_materials[f"{self.tcga}_protein_to_gene_dict"])

        # 3. Adding `NUM_ELASPIC_INTERFACE_ENTRIES` Column
        log.debug(f"Adding `NUM_ELASPIC_INTERFACE_ENTRIES` column ..")
        add_num_elaspic_interface_entries(
            preliminary_data, self.data_materials[f"{self.tcga}_protein_to_num_elaspic_interface_entries"]
        )

        # 4. Adding `NUM_DISRUPTIVE_ENTRIES` Column
        log.debug(f"Adding `NUM_DISRUPTIVE_ENTRIES` column ..")
        add_num_disruptive_entries(
            preliminary_data, self.data_materials[f"{self.tcga}_prediction_data"]
        )

        # 5. Adding `NUM_INCR_NOEFF_ENTRIES` Column
        log.debug(f"Adding `NUM_INCR_NOEFF_ENTRIES` column ..")
        add_num_incr_noeff_entries(
            preliminary_data, self.data_materials[f"{self.tcga}_prediction_data"]
        )

        # 6. Adding `NUM_UNIQUE_INTERACTORS` Column
        log.debug(f"Adding `NUM_UNIQUE_INTERACTORS` column ..")
        add_num_unique_interactors(
            preliminary_data, self.data_materials[f"{self.tcga}_protein_to_num_unique_interactors"]
        )

        # 7. Adding `PATIENT_CORE_COUNT` Column
        log.debug(f"Adding `PATIENT_CORE_COUNT` column ..")
        add_patient_core_count(
            preliminary_data,
            self.data_materials[f"{self.tcga}_snv_data_simplified"],
            elaspic_interface_data=self.data_materials[f"{self.tcga}_elaspic_interface_processed_data"],
            elaspic_core_data=self.data_materials[f"{self.tcga}_elaspic_core_data"]
        )

        # 8. Adding `PATIENT_INTERFACE_COUNT` Column
        log.debug(f"Adding `PATIENT_INTERFACE_COUNT` column ..")
        add_patient_interface_count(
            preliminary_data,
            self.data_materials[f"{self.tcga}_snv_data_simplified"],
            elaspic_interface_data=self.data_materials[f"{self.tcga}_elaspic_interface_processed_data"],
            elaspic_core_data=self.data_materials[f"{self.tcga}_elaspic_core_data"]
        )

        # update
        # 9. Adding `NUM_INTERFACE_PATIENTS_DISRUPTIVE_INTERACTOR` Column.
        log.debug(f"Adding `NUM_INTERFACE_PATIENTS_DISRUPTIVE_INTERACTOR` column ..")
        add_num_interface_patients_disruptive_interactor(
            preliminary_data,
            self.data_materials["patient_interaction_data"],
        )

        # 10. Adding `BASELINE` and `OUR_METHOD` Columns
        log.debug(f"Adding `BASELINE` and `OUR_METHOD` columns ..")
        proteins_to_counts_baseline_dict, proteins_to_counts_our_method_dict = counts_baseline_vs_our_method(
            proteins=self.data_materials[f"{self.tcga}_elaspic_proteins"],
            patients=self.data_materials[f"{self.tcga}_patients"],
            snv_data=self.data_materials[f"{self.tcga}_snv_data_simplified"],
            elaspic_core_data=self.data_materials[f"{self.tcga}_elaspic_core_data"],
            elaspic_interface_data=self.data_materials[f"{self.tcga}_elaspic_interface_processed_data"],
            prediction_data=self.data_materials[f"{self.tcga}_prediction_data"],
            add_core_flag_1_case_dict=None
        )
        add_baseline(preliminary_data, proteins_to_counts_baseline_dict)
        add_our_method(preliminary_data, proteins_to_counts_our_method_dict)

        # 11. Adding `OUR_METHOD / BASELINE` Column
        log.debug(f"Adding `OUR_METHOD / BASELINE` column ..")
        preliminary_data["OUR_METHOD/BASELINE"] = preliminary_data["OUR_METHOD"] / preliminary_data["BASELINE"]

        # 12. ELASPIC_COVERAGE
        log.debug(f"Adding `ELASPIC_COVERAGE` column ..")
        add_elaspic_coverage(
            preliminary_data,
            self.data_materials[f"{self.tcga}_elaspic_core_and_interface_data"],
            self.data_materials[f"{self.tcga}_snv_data_simplified"]
        )

        # 13. Adding Reference Dataset Columns: General and Cohort Specific
        # TODO Code refactor: add_reference_data_status
        log.debug(f"Adding Reference Dataset Columns: General and Cohort Specific columns ..")
        if self.reference_data_name == ReferenceDataset.CANCERMINE:
            add_cancermine_status(
                preliminary_data,
                cancermine_genes=self.data_materials["cancermine_all_genes"],
                cancermine_cohort_genes=self.data_materials[f"cancermine_{self.tcga}_genes"],
                tcga_type=self.tcga.upper()
            )

        elif self.reference_data_name == ReferenceDataset.CGC:
            add_cgc_status(
                preliminary_data,
                cgc_genes=self.data_materials["cgc_all_genes"],
                cgc_cohort_genes=self.data_materials[f"cgc_{self.tcga}_genes"],
                tcga_type=self.tcga.upper()
            )

        else:
            raise ValueError("Invalid Reference Dataset.")

        # Finally, assign the attribute.
        self.data_materials[f"{self.tcga}_preliminary_data"] = preliminary_data

        log.debug(f"{self.tcga}_preliminary_data is constructed.")

    def run_roc_curve_analysis(
            self,
            preliminary_data_name: str,
            state_variables=List[str]
    ):
        log.debug("Plotting ROC Curves ..")
        preliminary_data = self.data_materials[f"{preliminary_data_name}"].copy()

        ref_gene_column = state_variables[0]
        ref_gene_column_cohort = state_variables[1]

        self.auc_scores["default"] = roc_curve_analysis(
            reference_data_name=self.reference_data_name,
            preliminary_data=preliminary_data,
            ref_gene_column=ref_gene_column,
            cohort_specific=None
        )

        self.auc_scores[f"default_{self.tcga}"] = roc_curve_analysis(
            reference_data_name=self.reference_data_name,
            preliminary_data=preliminary_data,
            ref_gene_column=ref_gene_column_cohort,
            cohort_specific=self.tcga
        )

        # Baseline Non-Zero
        preliminary_data_baseline_nonzero = preliminary_data[preliminary_data['BASELINE'] != 0].copy()

        # Todo: pass another arg for prefix in generated table: e.g. "BASELINE_NON_ZERO"
        self.auc_scores[f"baseline_nonzero"] = roc_curve_analysis(
            reference_data_name=self.reference_data_name,
            preliminary_data=preliminary_data_baseline_nonzero,
            ref_gene_column=ref_gene_column,
            cohort_specific=None
        )

        self.auc_scores[f"baseline_nonzero_{self.tcga}"] = roc_curve_analysis(
            reference_data_name=self.reference_data_name,
            preliminary_data=preliminary_data_baseline_nonzero,
            ref_gene_column=ref_gene_column_cohort,
            cohort_specific=self.tcga
        )

    def export_auc_scores(self, file_name, overwrite=False):
        save_auc_scores(self.prediction_data_path, file_name, self.auc_scores, overwrite)

