from .conf import RECORD_FILE_PATH
from .log_script import ColorHandler

import os.path as op
import pandas as pd

import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')

log = logging.Logger("debug_runner", level=logging.DEBUG)
log.addHandler(ColorHandler())


def read_data(data_path):
    log.debug("Reading the record data ..")
    return pd.read_csv(data_path)


def is_file_exits(path):
    return op.isfile(path)


def save_record_data(data, data_path):
    log.debug("Saving record (written on disk).")
    data.to_csv(data_path, index=False)


class Entry:
    def __init__(
            self,
            gene=None,
            residue_position=None,
            submitted=None,
            downloaded=None,
            error_encountered=None,
            most_significant_codon_tier=None,
            CGC_status=None,
            url=None,
    ):
        self.gene = gene
        self.residue_position = residue_position
        self.error_encountered = error_encountered  # None
        self.submitted = submitted  # 1
        self.downloaded = downloaded  # 1
        self.most_significant_codon_tier = most_significant_codon_tier # TIER_2
        self.CGC_status = CGC_status  # Tier 1
        self.url = url

    def get_entry_data(self):
        entry_data = pd.DataFrame(
            {
                "GENE": [self.gene],
                "RESIDUE_POSITION": [self.residue_position],
                "ERROR_ENCOUNTERED": [self.error_encountered],
                "SUBMITTED": [self.submitted],
                "DOWNLOADED": [self.downloaded],
                "CGC_STATUS": [self.CGC_status],
                "MOST_SIGNIFICANT_CODON_TIER": [self.most_significant_codon_tier],
                "URL": [self.url],
            }
        )

        return entry_data


class Record:
    def __init__(self):
        self.record_data = None
        self.record_data_path = None

        self.load_record_data()

    def load_record_data(self):
        self.record_data_path = RECORD_FILE_PATH
        if is_file_exits(self.record_data_path):
            self.record_data = read_data(self.record_data_path)
        else:
            self.create_record_data()
            self.record_data = read_data(self.record_data_path)

    def create_record_data(self):
        log.info("Creating an empty record data ..")
        # Create a dummy data to get column names
        column_names = Entry().get_entry_data().columns
        record_data = pd.DataFrame(columns=column_names)

        save_record_data(record_data, self.record_data_path)

    def add_entry(self, entry: Entry):
        entry_data = entry.get_entry_data()
        self.record_data = self.record_data.append(entry_data, ignore_index=True)
        save_record_data(self.record_data, self.record_data_path)

    def is_already_recorded(self, gene, residue_position):
        search_data = self.record_data[
            (self.record_data["GENE"] == gene) &
            (self.record_data["RESIDUE_POSITION"] == int(residue_position))
        ]

        if search_data.empty:
            return False

        else:
            assert len(search_data) == 1
            return True
