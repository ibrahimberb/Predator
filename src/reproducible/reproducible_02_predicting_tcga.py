import papermill as pm
import os
import os.path as op

import sys
sys.path.insert(0, "../../")

FILE_NAME = op.basename(__file__)

# You may modify with it newly extracted model in reproducable_01_training_predator.py.
# e.g. PREDATOR_MODEL_PATH = "PredatorModels/PredatorModel_2022-06-24/e7935250/predator.pkl"
PREDATOR_MODEL_PATH = "PredatorModels/PredatorModel_2022-06-16/cc84a54e/predator.pkl"

if __name__ == "__main__":

    print(f"\nExecuting {FILE_NAME} ..\n")

    os.chdir("../")

    # Predicting on TCGA cohorts
    TCGA_LIST = ["BRCA", "COAD", "ESCA", "GBM", "HNSC", "OV"]
    for tcga in TCGA_LIST:
        print(f"Executing {tcga}")
        # Execute the PredatorStudy_<TCGA>
        pm.execute_notebook(
            f"PredatorStudy_{tcga}.ipynb",
            f"reproducible/Reproduced_PredatorStudy_{tcga}.ipynb",
            kernel_name="Predator", 
            parameters=dict(PREDATOR_MODEL_PATH=PREDATOR_MODEL_PATH),
            log_output=True,
            autosave_cell_every=60,
            stdout_file=sys.stdout,
            stderr_file=sys.stderr,
        )

    print(f"\n{FILE_NAME} is completed.\n")


 