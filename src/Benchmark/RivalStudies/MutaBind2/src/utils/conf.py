import pathlib

test_mode = False

MAX_RAM_ALLOWED_PERCENTAGE: int = 98

import os.path as op

COMPUTATION_TIME_ALLOWED = 1  # 3 # 10  # in seconds.

# Paths
DRIVER_PATH = r"C:\webdrivers\geckodriver.exe"
TEMP_DOWNLOAD_FOLDER_PATH = op.join("..", pathlib.Path().resolve(), "Firefox_download")
MUTABIND2_URL = "https://lilab.jysw.suda.edu.cn/research/mutabind2/"
PDB_FILES_PATH = op.join(pathlib.Path().resolve(), "..", "data", "pdb_files")
RECORD_PATH = "Record"
TEMP_FILES_PATH = op.join(pathlib.Path().resolve(), "temporary_files")

HEADLESS = False

# Firefox Log output disabled
PATH_TO_DEV_NULL = 'nul'
