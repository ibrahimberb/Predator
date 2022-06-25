import papermill as pm
import os
import os.path as op
import sys


sys.path.insert(0, "../../")

from src.helpers.helpers_predator.common import unzip_snv_files

FILE_NAME = op.basename(__file__)

# Feel free to change the parameters.
# Predictions

BRCA_PREDICTION_ID = "ed35a3a3/"
BRCA_PREDICTIONS_COMMON_PATH = "../data/predictions_datasets/brca_prediction_2022-06-17/" + BRCA_PREDICTION_ID
PREDICTION_BRCA_REDUCED_PATH = BRCA_PREDICTIONS_COMMON_PATH + "predictions_soft_2022-06-17.csv"

COAD_PREDICTION_ID = "84fd283a/"
COAD_PREDICTIONS_COMMON_PATH = "../data/predictions_datasets/coad_prediction_2022-06-17/" + COAD_PREDICTION_ID
PREDICTION_COAD_REDUCED_PATH = COAD_PREDICTIONS_COMMON_PATH + "predictions_soft_2022-06-17.csv"

ESCA_PREDICTION_ID = "f2d1e99a/"
ESCA_PREDICTIONS_COMMON_PATH = "../data/predictions_datasets/esca_prediction_2022-06-17/" + ESCA_PREDICTION_ID
PREDICTION_ESCA_REDUCED_PATH = ESCA_PREDICTIONS_COMMON_PATH + "predictions_soft_2022-06-17.csv"

GBM_PREDICTION_ID = "8d7f7caa/"
GBM_PREDICTIONS_COMMON_PATH = "../data/predictions_datasets/gbm_prediction_2022-06-17/" + GBM_PREDICTION_ID
PREDICTION_GBM_REDUCED_PATH = GBM_PREDICTIONS_COMMON_PATH + "predictions_soft_2022-06-17.csv"

HNSC_PREDICTION_ID = "76f498d9/"
HNSC_PREDICTIONS_COMMON_PATH = "../data/predictions_datasets/hnsc_prediction_2022-06-17/" + HNSC_PREDICTION_ID
PREDICTION_HNSC_REDUCED_PATH = HNSC_PREDICTIONS_COMMON_PATH + "predictions_soft_2022-06-17.csv"

OV_PREDICTION_ID = "865d1897/"
OV_PREDICTIONS_COMMON_PATH = "../data/predictions_datasets/ov_prediction_2022-06-17/" + OV_PREDICTION_ID
PREDICTION_OV_REDUCED_PATH = OV_PREDICTIONS_COMMON_PATH + "predictions_soft_2022-06-17.csv"

# ELASPIC Core and Interface datasets
ELASPIC_RESULTS_COMMON_PATH = "../data/Elaspic_merged_results/"  # elaspic_results_datasets

BRCA_CORE_PATH = ELASPIC_RESULTS_COMMON_PATH + "BRCA_Core_2021-11-17.txt"
BRCA_INTERFACE_PATH = ELASPIC_RESULTS_COMMON_PATH + "BRCA_Interface_2021-11-17.txt"

COAD_CORE_PATH = ELASPIC_RESULTS_COMMON_PATH + "COAD_Core_2022-01-06.txt"
COAD_INTERFACE_PATH = ELASPIC_RESULTS_COMMON_PATH + "COAD_Interface_2022-01-06.txt"

ESCA_CORE_PATH = ELASPIC_RESULTS_COMMON_PATH + "ESCA_Core_2021-11-17.txt"
ESCA_INTERFACE_PATH = ELASPIC_RESULTS_COMMON_PATH + "ESCA_Interface_2021-11-17.txt"

GBM_CORE_PATH = ELASPIC_RESULTS_COMMON_PATH + "GBM_Core_2021-11-17.txt"
GBM_INTERFACE_PATH = ELASPIC_RESULTS_COMMON_PATH + "GBM_Interface_2021-11-17.txt"

HNSC_CORE_PATH = ELASPIC_RESULTS_COMMON_PATH + "HNSC_Core_2021-11-17.txt"
HNSC_INTERFACE_PATH = ELASPIC_RESULTS_COMMON_PATH + "HNSC_Interface_2021-11-17.txt"

OV_CORE_PATH = ELASPIC_RESULTS_COMMON_PATH + "OV_Core_2021-11-17.txt"
OV_INTERFACE_PATH = ELASPIC_RESULTS_COMMON_PATH + "OV_Interface_2021-11-17.txt"


# SNV_PATHS
SNV_COMMON_PATH = "../data/snv_datasets"
BRCA_SNV_PATH = op.join(SNV_COMMON_PATH, "SNV_BRCA_hg38_2021-09-22.csv")
COAD_SNV_PATH = op.join(SNV_COMMON_PATH, "SNV_COAD_hg38_2021-09-22.csv")
ESCA_SNV_PATH = op.join(SNV_COMMON_PATH, "SNV_ESCA_hg38_2021-09-22.csv")
GBM_SNV_PATH = op.join(SNV_COMMON_PATH, "SNV_GBM_hg38_2021-09-22.csv")
HNSC_SNV_PATH = op.join(SNV_COMMON_PATH, "SNV_HNSC_hg38_2021-09-22.csv")
OV_SNV_PATH = op.join(SNV_COMMON_PATH, "SNV_OV_hg38_2021-09-22.csv")

PARAMETERS = dict(
    PREDICTION_BRCA_REDUCED_PATH=PREDICTION_BRCA_REDUCED_PATH,
    PREDICTION_COAD_REDUCED_PATH=PREDICTION_COAD_REDUCED_PATH,
    PREDICTION_ESCA_REDUCED_PATH=PREDICTION_ESCA_REDUCED_PATH,
    PREDICTION_GBM_REDUCED_PATH=PREDICTION_GBM_REDUCED_PATH,
    PREDICTION_HNSC_REDUCED_PATH=PREDICTION_HNSC_REDUCED_PATH,
    PREDICTION_OV_REDUCED_PATH=PREDICTION_OV_REDUCED_PATH,
    BRCA_CORE_PATH=BRCA_CORE_PATH,
    BRCA_INTERFACE_PATH=BRCA_INTERFACE_PATH,
    COAD_CORE_PATH=COAD_CORE_PATH,
    ESCA_CORE_PATH=ESCA_CORE_PATH,
    ESCA_INTERFACE_PATH=ESCA_INTERFACE_PATH,
    GBM_CORE_PATH=GBM_CORE_PATH,
    GBM_INTERFACE_PATH=GBM_INTERFACE_PATH,
    HNSC_CORE_PATH=HNSC_CORE_PATH,
    HNSC_INTERFACE_PATH=HNSC_INTERFACE_PATH,
    OV_CORE_PATH=OV_CORE_PATH,
    OV_INTERFACE_PATH=OV_INTERFACE_PATH,
    SNV_COMMON_PATH=SNV_COMMON_PATH,
    BRCA_SNV_PATH=BRCA_SNV_PATH,
    COAD_SNV_PATH=COAD_SNV_PATH,
    ESCA_SNV_PATH=ESCA_SNV_PATH,
    GBM_SNV_PATH=GBM_SNV_PATH,
    HNSC_SNV_PATH=HNSC_SNV_PATH,
    OV_SNV_PATH=OV_SNV_PATH,
)

if __name__ == "__main__":
    
    print(f"\nExecuting {FILE_NAME} ..\n")

    unzip_snv_files("../../data/snv_datasets/")

    os.chdir("../")

    pm.execute_notebook(
        f"Disruptive_interactions_per_patient.ipynb",
        f"reproducible/Reproduced_Disruptive_interactions_per_patient.ipynb",
        kernel_name="Predator", 
        # parameters=PARAMETERS,
        log_output=True,
        autosave_cell_every=60,
        stdout_file=sys.stdout,
        stderr_file=sys.stderr,
    )

    print(f"\n{FILE_NAME} is completed.\n")


 