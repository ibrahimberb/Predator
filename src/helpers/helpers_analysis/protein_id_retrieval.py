from pandas import read_csv
import os.path as op

# for PyCharm
# from src.helpers.mylogger import get_handler

from ..mylogger import get_handler
import logging

handler = get_handler(log_type="module")

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.WARNING)


class ProteinIDFetcher:
    """Fetch only. Don't download."""
    DATA_DIR = "helpers/helpers_analysis/gene_retrieval"
    # MAPPING_DATA_NAME = op.join(DATA_DIR, "UNIPROT_GENE_MAPPING.csv")

    def __init__(self, mapping_data_path=None):
        if mapping_data_path is None:
            self.MAPPING_DATA_NAME = op.join(self.DATA_DIR, "UNIPROT_GENE_MAPPING.csv")
        else:
            self.MAPPING_DATA_NAME = mapping_data_path

        self.mapping_data = None
        self.load_mapping_data()

    def load_mapping_data(self):
        if op.isfile(self.MAPPING_DATA_NAME):
            log.info("Reading mapping data..")
            mapping_data = read_csv(self.MAPPING_DATA_NAME, index_col="UNIPROT")
        else:
            log.critical("No mapping data found. Make sure you have the mapping data.")
            raise FileNotFoundError("No mapping data found!")

        self.mapping_data = mapping_data
        log.info("Mapping data loaded.")

    def fetch(self, gene):
        log.debug(f"Fetching gene id of protein {gene}")

        try:
            protein = self.mapping_data[
                self.mapping_data["GENE"] == gene
            ].index.to_list()

        except ValueError:
            raise

        return protein
