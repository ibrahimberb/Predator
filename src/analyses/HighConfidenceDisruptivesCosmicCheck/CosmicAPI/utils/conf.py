import pathlib
import os.path as op

test_mode = False

ALLOWED_RAM_PERCENTAGE: int = 93

# Paths
DRIVER_PATH = r"C:\webdrivers\geckodriver.exe"
TEMP_DOWNLOAD_FOLDER_PATH = op.join("..", pathlib.Path().resolve(), "Firefox_download")
COSMIC_URL = "http://www.google.com"
PDB_FILES_PATH = op.join(pathlib.Path().resolve(), "..", "data", "pdb_files")
RECORD_FILE_NAME = "record_data.csv"
RECORD_FOLDER_PATH = "Record"
RECORD_FILE_PATH = op.join(RECORD_FOLDER_PATH, RECORD_FILE_NAME)
DATA_FOLDER_PATH = "Data"

USER_CREDENTIALS_PATH = "account_credentials.txt"

HEADLESS = True

# Firefox Log output disabled
PATH_TO_DEV_NULL = 'nul'
