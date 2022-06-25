import papermill as pm
import os
import os.path as op

import sys
sys.path.insert(0, "../../")

FILE_NAME = op.basename(__file__)

# You may provide your own parameters and give as `parameters` dictionary below.
# Predictions paths
PREDICTION_BRCA_REDUCED_PATH = "../data/predictions_datasets/brca_prediction_2022-06-24/4f29e8e7/predictions_soft_2022-06-24.csv"
PREDICTION_COAD_REDUCED_PATH = "../data/predictions_datasets/coad_prediction_2022-06-24/b11e0844/predictions_soft_2022-06-24.csv"
PREDICTION_ESCA_REDUCED_PATH = "../data/predictions_datasets/esca_prediction_2022-06-24/ee706abc/predictions_soft_2022-06-24.csv"
PREDICTION_GBM_REDUCED_PATH = "../data/predictions_datasets/gbm_prediction_2022-06-24/9917e417/predictions_soft_2022-06-24.csv"
PREDICTION_HNSC_REDUCED_PATH = "../data/predictions_datasets/hnsc_prediction_2022-06-24/00c91c32/predictions_soft_2022-06-24.csv"
PREDICTION_OV_REDUCED_PATH = "../data/predictions_datasets/ov_prediction_2022-06-24/5bec840a/predictions_soft_2022-06-24.csv"

# Patient interaction data paths
BRCA_PATIENT_INTERACTION_DATA_PATH = "../data/patient_interaction_datasets/BRCA_patient_interactions_analysis_table_2022-06-24.xlsx"
COAD_PATIENT_INTERACTION_DATA_PATH = "../data/patient_interaction_datasets/COAD_patient_interactions_analysis_table_2022-06-24.xlsx"
ESCA_PATIENT_INTERACTION_DATA_PATH = "../data/patient_interaction_datasets/ESCA_patient_interactions_analysis_table_2022-06-24.xlsx"
GBM_PATIENT_INTERACTION_DATA_PATH = "../data/patient_interaction_datasets/GBM_patient_interactions_analysis_table_2022-06-24.xlsx"
HNSC_PATIENT_INTERACTION_DATA_PATH = "../data/patient_interaction_datasets/HNSC_patient_interactions_analysis_table_2022-06-24.xlsx"
OV_PATIENT_INTERACTION_DATA_PATH = "../data/patient_interaction_datasets/OV_patient_interactions_analysis_table_2022-06-24.xlsx"

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
    BRCA_PATIENT_INTERACTION_DATA_PATH=BRCA_PATIENT_INTERACTION_DATA_PATH,
    COAD_PATIENT_INTERACTION_DATA_PATH=COAD_PATIENT_INTERACTION_DATA_PATH,
    ESCA_PATIENT_INTERACTION_DATA_PATH=ESCA_PATIENT_INTERACTION_DATA_PATH,
    GBM_PATIENT_INTERACTION_DATA_PATH=GBM_PATIENT_INTERACTION_DATA_PATH,
    HNSC_PATIENT_INTERACTION_DATA_PATH=HNSC_PATIENT_INTERACTION_DATA_PATH,
    OV_PATIENT_INTERACTION_DATA_PATH=OV_PATIENT_INTERACTION_DATA_PATH,
    BRCA_SNV_PATH=BRCA_SNV_PATH,
    COAD_SNV_PATH=COAD_SNV_PATH,
    ESCA_SNV_PATH=ESCA_SNV_PATH,
    GBM_SNV_PATH=GBM_SNV_PATH,
    HNSC_SNV_PATH=HNSC_SNV_PATH,
    OV_SNV_PATH=OV_SNV_PATH,
)

if __name__ == "__main__":

    print(f"\nExecuting {FILE_NAME} ..\n")

    os.chdir("../")

    # Predicting on TCGA cohorts
    TCGA_LIST = ["BRCA", "COAD", "ESCA", "GBM", "HNSC", "OV"]
    for tcga in TCGA_LIST:
        print(f"Executing {tcga}")
        # Execute the PredatorAnalysis_<TCGA>_CGC
        pm.execute_notebook(
            f"PredatorAnalysis_{tcga}_CGC.ipynb",
            f"reproducible/Reproduced_PredatorAnalysis_{tcga}.ipynb",
            kernel_name="Predator", 
            # parameters=PARAMETERS,
            log_output=True,
            autosave_cell_every=60,
            stdout_file=sys.stdout,
            stderr_file=sys.stderr,
        )

    print(f"\n{FILE_NAME} is completed.\n")


 