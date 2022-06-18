from pathlib import Path
import pandas as pd

from .preprocessing import (
    process_snv,
    shorten_patient_ids,
    simplify_snv_data,
    simplify_elaspic_data,
    remove_duplicated_entries,
    convert_primary_isomer,
)

from ..mylogger import get_handler
import logging

handler = get_handler('module')

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class ReferenceDataset:
    CANCERMINE = 'cancermine'
    CGC = 'cgc'


def load_snv_datasets(tcga: str, snv_path: Path, data_materials: dict):
    log.info(f"Loading {tcga} SNV datasets ..")
    tcga_snv_data = pd.read_csv(snv_path, low_memory=False)
    log.debug(f"{tcga} SNV data size: {tcga_snv_data.shape}")

    tcga_snv_data_processed = process_snv(tcga_snv_data)
    log.debug(f"{tcga} SNV data processed size: {tcga_snv_data_processed.shape}")

    tcga_snv_data_processed = shorten_patient_ids(tcga_snv_data_processed)
    tcga_snv_data_simplified = simplify_snv_data(tcga_snv_data_processed)

    data_materials[f"{tcga}_snv_data"] = tcga_snv_data
    data_materials[f"{tcga}_snv_data_processed"] = tcga_snv_data_processed
    data_materials[f"{tcga}_snv_data_simplified"] = tcga_snv_data_simplified

    log.info(f"{tcga} SNV datasets are loaded.")


def load_elaspic_core_data(tcga, elaspic_core_path, data_materials):
    # Core
    log.debug("Loading ELASPIC CORE data materials ..")

    tcga_elaspic_core_data = pd.read_csv(elaspic_core_path, sep='\t', low_memory=False)
    tcga_elaspic_core_data.drop_duplicates(keep="first", inplace=True)

    log.debug(f"{tcga} ELASPIC CORE data size: {tcga_elaspic_core_data.shape}")

    tcga_elaspic_core_data_simplified = simplify_elaspic_data(tcga_elaspic_core_data)
    tcga_elaspic_core_data_simplified.drop_duplicates(keep="first", inplace=True)

    assert tcga_elaspic_core_data_simplified.duplicated().sum() == 0  # No duplicated entries in ELASPIC core data.

    data_materials[f"{tcga}_elaspic_core_data"] = tcga_elaspic_core_data
    data_materials[f"{tcga}_elaspic_core_data_simplified"] = tcga_elaspic_core_data_simplified


def load_elaspic_interface_data(tcga, elaspic_interface_path, data_materials):
    # Interface
    log.debug("Loading ELASPIC INTERFACE data materials ..")

    tcga_elaspic_interface_data = pd.read_csv(elaspic_interface_path, sep='\t', low_memory=False)
    tcga_elaspic_interface_data.drop_duplicates(keep="first", inplace=True)

    log.debug(f"{tcga} ELASPIC INTERFACE data size: {tcga_elaspic_interface_data.shape}")

    tcga_elaspic_interface_processed_data = simplify_elaspic_data(tcga_elaspic_interface_data)
    tcga_elaspic_interface_processed_data = convert_primary_isomer(
        "Interactor_UniProt_ID", tcga_elaspic_interface_processed_data
    )
    tcga_elaspic_interface_processed_data = remove_duplicated_entries(tcga_elaspic_interface_processed_data)
    tcga_elaspic_interface_processed_data.reset_index(drop=True, inplace=True)

    data_materials[f"{tcga}_elaspic_interface_data"] = tcga_elaspic_interface_data
    data_materials[f"{tcga}_elaspic_interface_processed_data"] = tcga_elaspic_interface_processed_data


def load_elaspic_core_and_interface_combined_data(tcga, core_data, interface_data, data_materials):
    # Core and Interface
    log.debug(f"Loading {tcga} ELASPIC CORE & INTERFACE COMBINED data ..")

    log.debug(f"\nELASPIC CORE DATA"
              f"\n{core_data.shape}"
              f"\n{core_data.head(3)}")

    log.debug(f'\nELASPIC INTERFACE DATA'
              f'\n{interface_data.shape}'
              f'\n{interface_data.head(3)}')

    tcga_elaspic_core_and_interface_data = pd.concat([core_data, interface_data])

    log.debug(f"Combined data shape: {tcga_elaspic_core_and_interface_data.shape}")

    data_materials[f"{tcga}_elaspic_core_and_interface_data"] = tcga_elaspic_core_and_interface_data


def load_elaspic_datasets(tcga, elaspic_core_path, elaspic_interface_path, data_materials):
    log.info(f"Loading {tcga} ELASPIC datasets ..")

    load_elaspic_core_data(tcga, elaspic_core_path, data_materials)
    load_elaspic_interface_data(tcga, elaspic_interface_path, data_materials)
    load_elaspic_core_and_interface_combined_data(
        tcga,
        data_materials[f"{tcga}_elaspic_core_data_simplified"],
        data_materials[f"{tcga}_elaspic_interface_processed_data"],
        data_materials
    )

    log.info(f"{tcga} ELASPIC datasets are loaded.")


def load_prediction_dataset(tcga, prediction_data_path, data_materials):
    log.info(f"Loading {tcga} Prediction dataset ..")

    # Prediction reduced already
    tcga_prediction_data = pd.read_csv(prediction_data_path, low_memory=False)
    log.debug(f"{tcga} Prediction data shape: {tcga_prediction_data.shape}")

    data_materials[f"{tcga}_prediction_data"] = tcga_prediction_data

    log.info(f"{tcga} Prediction dataset is loaded.")


def load_reference_dataset(
        tcga,
        reference_data_name,
        reference_data_spec_cohort_path,
        reference_data_path,
        data_materials
):
    log.info(f"Loading reference dataset {reference_data_name} for TCGA {tcga} ..")
    if reference_data_name == ReferenceDataset.CANCERMINE:
        load_cancermine_genes_cohort(reference_data_spec_cohort_path, data_materials, cohort=tcga)
        load_cancermine_genes_cohort(reference_data_path, data_materials, cohort='all')

        log.info(f"{reference_data_name} reference datasets loaded.")

    elif reference_data_name == ReferenceDataset.CGC:
        load_cgc_genes_cohort(reference_data_spec_cohort_path, data_materials, cohort=tcga)
        load_cgc_genes_cohort(reference_data_path, data_materials, cohort='all')

    else:
        raise ValueError(f"Invalid name `{reference_data_name}`")


def load_cancermine_genes_cohort(reference_data_path, data_materials, cohort):
    log.debug(f"Loading Cancermine reference dataset cohort {cohort} ..")

    with open(reference_data_path, 'r') as genes_file:
        cancermine_genes = [line.strip() for line in genes_file.readlines()]

    log.debug(f"Number of genes: {len(cancermine_genes)}")
    log.debug(f"First five genes: {cancermine_genes[:5]}")

    data_materials[f"cancermine_{cohort}_genes"] = cancermine_genes


# fixme: refactor (duplicated code: genes txt files are in the same format, so we can use the same function.)
def load_cgc_genes_cohort(reference_data_path, data_materials, cohort):
    log.debug(f"Loading CGC reference dataset cohort {cohort} ..")

    with open(reference_data_path, 'r') as genes_file:
        cgc_genes = [line.strip() for line in genes_file.readlines()]

    log.debug(f"Number of genes: {len(cgc_genes)}")
    log.debug(f"First five genes: {cgc_genes[:5]}")

    data_materials[f"cgc_{cohort}_genes"] = cgc_genes


def load_patient_interaction_dataset(patient_interaction_data_path, data_materials):
    log.debug(f"Loading Patient Interaction Dataset ..")
    patient_interaction_data = pd.read_excel(patient_interaction_data_path)
    data_materials["patient_interaction_data"] = patient_interaction_data
