from .conf import RECORD_PATH
from .log_script import ColorHandler

import os.path as op
import pandas as pd

import logging

from .misc import amino_acid_shorten

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


def post_check_entry_data(data):
    # Check there's only one entry
    assert len(data) == 1

    # Check if mutation is correct.
    [mutation] = data["MUTATION"]
    [mutation_registered] = data["MUTATION_REGISTERED"]
    assert mutation == mutation_registered

    # check if chain is correct.
    [chain] = data["CHAIN_ID_1"]
    [post_chain] = data["MUTATED_CHAIN"]
    assert post_chain == chain


class Entry:
    def __init__(
            self,
            mutation_effect_label=None,
            protein=None,
            mutation=None,
            interactor=None,
            pdb_id=None,
            chain_id_1=None,
            chain_id_2=None,
            submitted=None,
            saved=None,
            error_encountered=None,
            result_url=None,
            mutation_details=None,
    ):
        self.mutation_effect_label = mutation_effect_label
        self.protein = protein  # "Q9BPZ3"
        self.mutation = mutation  # "F118A"
        self.interactor = interactor  # "P11940"
        self.pdb_id = pdb_id  # "1jgn"
        self.chain_id_1 = chain_id_1  # "B"
        self.chain_id_2 = chain_id_2  # "B"
        self.error_encountered = error_encountered  # None
        self.submitted = submitted  # 1
        self.saved = saved  # 1
        self.result_url = result_url  # "http://biosig.unimelb.edu.au/mcsm_ppi2/results_prediction/164486655328"

        # Mutation Details
        if mutation_details is None:
            self.job_id = None
            self.pdb_id_registered = None
            self.mutated_chain = None
            self.mutation_registered = None
            self.ddG_predicted = None
            self.interface = None
            self.deleterious = None

        else:
            self.job_id = mutation_details["job_id"]
            self.pdb_id_registered = mutation_details["pdb_id_registered"]
            self.mutated_chain = mutation_details["mutated_chain"]
            self.mutation_registered = mutation_details["mutation_registered"]
            self.ddG_predicted = mutation_details["ddG_predicted"]
            self.interface = mutation_details["interface"]
            self.deleterious = mutation_details["deleterious"]

    def get_entry_data(self):
        entry_data = pd.DataFrame(
            {
                "MUTATION_EFFECT_LABEL": [self.mutation_effect_label],
                "PROTEIN": [self.protein],
                "MUTATION": [self.mutation],
                "INTERACTOR": [self.interactor],
                "PDB_ID": [self.pdb_id],
                "CHAIN_ID_1": [self.chain_id_1],
                "CHAIN_ID_2": [self.chain_id_2],
                "ERROR_ENCOUNTERED": [self.error_encountered],
                "SUBMITTED": [self.submitted],
                "SAVED": [self.saved],
                "RESULT_URL": [self.result_url],
                "JOB_ID": [self.job_id],
                "PDB_ID_REGISTERED": [self.pdb_id_registered],
                "MUTATED_CHAIN": [self.mutated_chain],
                "MUTATION_REGISTERED": [self.mutation_registered],
                "DDG_PREDICTED": [self.ddG_predicted],
                "INTERFACE": [self.interface],
                "DELETERIOUS": [self.deleterious],
            }
        )

        return entry_data


class Record:
    def __init__(self):
        self.record_data = None
        self.record_data_path = None

        self.load_record_data()

    def load_record_data(self):
        record_data_path = op.join(RECORD_PATH, "record_data.csv")
        self.record_data_path = record_data_path
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

    def add_entry(self, entry: Entry, successful_entry):
        # merge two dataframes
        entry_data = entry.get_entry_data()
        if successful_entry:
            post_check_entry_data(entry_data)
        self.record_data = self.record_data.append(entry_data, ignore_index=True)
        save_record_data(self.record_data, self.record_data_path)


class TrainDataIndexes:
    def __init__(self):
        self.successful_train_indexes_file_path = op.join(RECORD_PATH, "successful_train_index.txt")
        self.failed_train_indexes_file_path = op.join(RECORD_PATH, "failed_train_index.txt")

        if is_file_exits(self.successful_train_indexes_file_path):
            self.successful_indexes = self.read_indexes(self.successful_train_indexes_file_path)
        else:
            self.create_index_file(self.successful_train_indexes_file_path)
            self.successful_indexes = self.read_indexes(self.successful_train_indexes_file_path)

        if is_file_exits(self.failed_train_indexes_file_path):
            self.failed_indexes = self.read_indexes(self.failed_train_indexes_file_path)
        else:
            self.create_index_file(self.failed_train_indexes_file_path)
            self.failed_indexes = self.read_indexes(self.failed_train_indexes_file_path)

    def add_successful_index(self, index):
        self.successful_indexes.append(index)
        self.write_successful_index_file(self.successful_train_indexes_file_path)

    def add_failed_index(self, index):
        self.failed_indexes.append(index)
        self.write_failed_index_file(self.failed_train_indexes_file_path)

    def write_successful_index_file(self, file_path):
        log.debug("Writing successful index file")
        with open(file_path, "w") as fout:
            for index in sorted(self.successful_indexes):
                fout.write(f"{index}\n")

    def write_failed_index_file(self, file_path):
        log.debug("Writing failed index file")
        with open(file_path, "w") as fout:
            for index in sorted(self.failed_indexes):
                fout.write(f"{index}\n")

    @staticmethod
    def create_index_file(file_path):
        log.info("Creating index file")
        with open(file_path, "w"):
            pass

    @staticmethod
    def read_indexes(file_path):
        with open(file_path, "r") as fin:
            lines = fin.readlines()
            indexes = [int(line) for line in lines]

        return indexes
