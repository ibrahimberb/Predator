import pathlib

test_mode = False

ALLOWED_RAM_PERCENTAGE: int = 93

import os.path as op

COMPUTATION_TIME_ALLOWED = 1  # 3 # 10  # in seconds.

# Paths
DRIVER_PATH = r"C:\webdrivers\geckodriver.exe"
TEMP_DOWNLOAD_FOLDER_PATH = op.join("..", pathlib.Path().resolve(), "Firefox_download")
MCSM_PPI_MANY_URL = "http://biosig.unimelb.edu.au/mcsm_ppi2/submit_prediction"
PDB_FILES_PATH = op.join(pathlib.Path().resolve(), "..", "data", "pdb_files")
RECORD_PATH = "Record"

HEADLESS = False

# Firefox Log output disabled
PATH_TO_DEV_NULL = 'nul'
