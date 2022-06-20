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
    [wildtype] = data["MUTATION_DETAILS_WILDTYPE"]
    wildtype = amino_acid_shorten(wildtype)
    [position] = data["MUTATION_DETAILS_POSITION"]
    [mutant] = data["MUTATION_DETAILS_MUTANT"]
    mutant = amino_acid_shorten(mutant)
    post_mutation = f"{wildtype}{position}{mutant}"
    assert mutation == post_mutation

    # check if chain is correct.
    [chain] = data["CHAIN_ID"]
    [post_chain] = data["MUTATION_DETAILS_CHAIN"]
    assert post_chain == chain

    # assert its change is correct.
    [change_value] = data["PREDICTED_AFFINITY_CHANGE_VALUE"]
    change_value = float(change_value)
    [change] = data["PREDICTED_AFFINITY_CHANGE"]

    if change_value < 0:
        change_converted = "Decreasing affinity"
    elif change_value > 0:
        change_converted = "Increasing affinity"
    else:
        log.critical(f"change_converted: {change_value}")
        # It seems that change value could be 0 also.
        raise

    assert change == change_converted


class Entry:
    def __init__(
            self,
            mutation_effect_label=None,
            protein=None,
            mutation=None,
            interactor=None,
            pdb_id=None,
            chain_id=None,
            submitted=None,
            saved=None,
            error_encountered=None,
            result_url=None,
            predicted_affinity_change_value=None,
            predicted_affinity_change=None,
            mutation_details=None,
    ):
        self.mutation_effect_label = mutation_effect_label
        self.protein = protein  # "Q9BPZ3"
        self.mutation = mutation  # "F118A"
        self.interactor = interactor  # "P11940"
        self.pdb_id = pdb_id  # "1jgn"
        self.chain_id = chain_id  # "B"
        self.error_encountered = error_encountered  # None
        self.submitted = submitted  # 1
        self.saved = saved  # 1
        self.result_url = result_url  # "http://biosig.unimelb.edu.au/mcsm_ppi2/results_prediction/164486655328"

        # Predicted result
        self.predicted_affinity_change_value = predicted_affinity_change_value  # -2.18
        self.predicted_affinity_change = predicted_affinity_change  # "Decreasing affinity"

        # Mutation Details
        if mutation_details is None:
            self.mutation_details_chain = None
            self.mutation_details_position = None
            self.mutation_details_wildtype = None
            self.mutation_details_mutant = None
            self.mutation_details_distance_from_closest_partner = None
        else:
            self.mutation_details_chain = mutation_details["chain"]  # "A"
            self.mutation_details_position = mutation_details["position"]  # 40
            self.mutation_details_wildtype = mutation_details["wildtype"]  # "TYR"
            self.mutation_details_mutant = mutation_details["mutant"]  # "CYS"
            self.mutation_details_distance_from_closest_partner = mutation_details["dis_closest_partner"]  # "3.494Å"

    def get_entry_data(self):
        entry_data = pd.DataFrame(
            {
                "MUTATION_EFFECT_LABEL": [self.mutation_effect_label],
                "PROTEIN": [self.protein],
                "MUTATION": [self.mutation],
                "INTERACTOR": [self.interactor],
                "PDB_ID": [self.pdb_id],
                "CHAIN_ID": [self.chain_id],
                "ERROR_ENCOUNTERED": [self.error_encountered],
                "SUBMITTED": [self.submitted],
                "SAVED": [self.saved],
                "RESULT_URL": [self.result_url],
                "PREDICTED_AFFINITY_CHANGE_VALUE": [self.predicted_affinity_change_value],
                "PREDICTED_AFFINITY_CHANGE": [self.predicted_affinity_change],
                "MUTATION_DETAILS_CHAIN": [self.mutation_details_chain],
                "MUTATION_DETAILS_POSITION": [self.mutation_details_position],
                "MUTATION_DETAILS_WILDTYPE": [self.mutation_details_wildtype],
                "MUTATION_DETAILS_MUTANT": [self.mutation_details_mutant],
                "MUTATION_DETAILS_DISTANCE_FROM_CLOSEST_PARTNER": [self.mutation_details_distance_from_closest_partner],
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
