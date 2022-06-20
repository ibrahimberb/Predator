from typing import Union
import time

from tqdm import tqdm
import os.path as op

from .conf import PDB_FILES_PATH

from .log_script import ColorHandler
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')

log = logging.Logger("debug_runner", level=logging.DEBUG)
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


def get_PDB_file_path(pdb_id):
    return op.join(PDB_FILES_PATH, f"{pdb_id}.pdb")


def amino_acid_shorten(aa_three_letter):
    if len(aa_three_letter) != 3:
        log.critical(f"aa_three_letter: {aa_three_letter}")
        raise ValueError('Input length should be a multiple of three')

    mapping = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
               'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
               'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
               'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    return mapping[aa_three_letter]
