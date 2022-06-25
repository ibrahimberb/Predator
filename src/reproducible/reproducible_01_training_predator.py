import papermill as pm
import os
import warnings

import os.path as op

from sklearn.exceptions import UndefinedMetricWarning

import sys
sys.path.insert(0, "../../")

warnings.filterwarnings(
    'ignore', category=RuntimeWarning, append=True
)

warnings.filterwarnings(
    'ignore', category=UndefinedMetricWarning, append=True
)

FILE_NAME = op.basename(__file__)

N_CORES = 4

if __name__ == "__main__":

    print(f"\nExecuting {FILE_NAME} ..\n")

    os.chdir("../")

    # Execute the PredatorStudyModel
    pm.execute_notebook(
        "PredatorStudyModel.ipynb",
        "reproducible/Reproduced_PredatorStudyModel.ipynb",
        kernel_name="Predator", 
        parameters=dict(N_CORES=N_CORES),
        log_output=True,
        autosave_cell_every=60,
        stdout_file=sys.stdout,
        stderr_file=sys.stderr,
    )

    print(f"\n{FILE_NAME} is completed.\n")


 