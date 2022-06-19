import re
from datetime import datetime
from pathlib import Path
from typing import Union
import time

from tqdm import tqdm

import os.path as op

from .conf import DATA_FOLDER_PATH
from .log_script import ColorHandler, MyLog
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')

log = MyLog("debug_runner", level=logging.DEBUG)
log.addHandler(ColorHandler())

Seconds = Union[int, float]


def wait(duration: Seconds, desc=None):
    if desc is None:
        desc = "[WAIT_DELAY]"
    if duration == 0:
        return
    if duration <= 1:
        time.sleep(duration)

    elif isinstance(duration, int):
        for _ in tqdm(range(duration), desc=desc, position=0, leave=True):
            time.sleep(1)

    else:
        time.sleep(duration)


def load_credentials(credentials_path):
    with open(credentials_path) as fin:
        lines = [line for line in fin.readlines() if not line.startswith("#")]

    email, password = lines

    credentials = {
        "email": email,
        "password": password,
    }

    return credentials


def get_cosmic_url(gene, mutation):
    residue_position = get_residue_position(mutation)
    url = fr"https://cancer.sanger.ac.uk/cmc/gene/{gene}/codon/{residue_position}"
    return url


def get_residue_position(mutation):
    aa_original, residue_position, aa_mutated = re.match(r"^([A-Z])(\d+)([A-Z])$", mutation).groups()
    return residue_position


def save_results(data, folder_path, file_name):
    current_date = datetime.now().strftime("%Y-%m-%d")
    Path(f"{folder_path}").mkdir(parents=True, exist_ok=True)
    file_name = op.join(folder_path, f"{file_name}_{current_date}.csv")
    data.to_csv(file_name, index=False)
    log.success(f"Data {file_name} exported successfully.")


def save_CGC_data(data, gene, residue_position):
    folder_path = op.join(DATA_FOLDER_PATH, gene)
    file_name = f"{gene}_{residue_position}"
    save_results(data=data, folder_path=folder_path, file_name=file_name)
