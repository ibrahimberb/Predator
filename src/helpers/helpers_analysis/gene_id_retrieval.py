import re
import time

import requests

from pandas import DataFrame
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


class GeneIDRetriever:
    DATA_DIR = "helpers/helpers_analysis/gene_retrieval"
    MAPPING_DATA_NAME = op.join(DATA_DIR, "UNIPROT_GENE_MAPPING.csv")

    def __init__(self):
        self.mapping_data = None
        self.load_mapping_data()

    def load_mapping_data(self):
        if op.isfile(self.MAPPING_DATA_NAME):
            log.info("Reading mapping data..")
            mapping_data = read_csv(self.MAPPING_DATA_NAME, index_col="UNIPROT")
        else:
            log.info("Creating new mapping data..")
            mapping_data = DataFrame(columns=["UNIPROT", "GENE"])
            mapping_data.set_index('UNIPROT', inplace=True)
            # self.save_mapping_data()  <- this should belong here.

        self.mapping_data = mapping_data
        log.info("Mapping data loaded.")

    def save_mapping_data(self):
        log.info("Saving mapping data..")
        if self.mapping_data is None:
            raise ValueError("Mapping data is None.")
        self.mapping_data.to_csv(self.MAPPING_DATA_NAME)

    def fetch(self, protein):
        log.debug(f"Fetching gene id of protein {protein}")

        try:
            gene = self.mapping_data.loc[protein]["GENE"]

        except KeyError:
            log.info(f"Protein {protein} not found in database. Retrieving ..")
            gene = self._retrieve(protein)
            self.insert_entry(protein=protein, gene=gene)
            self.save_mapping_data()
            assert self.mapping_data.loc[protein]["GENE"] == gene  # Ensuring insertion is reflected.

        return gene

    def insert_entry(self, protein, gene):
        """
        Adds an entry containing protein and corresponding gene.
        """
        mapping_data = self.mapping_data.copy()
        assert protein not in mapping_data.index, f"Entry with protein {protein} already exists."
        # Adding entry that maps protein to its gene.
        mapping_data.loc[protein] = gene
        self.mapping_data = mapping_data
        log.debug(f"({protein}, {gene}) pairs are added.")

    def _retrieve(self, protein):
        """
        Retrieves the Gene name of given protein from UniProt API.
        """
        address = "http://www.uniprot.org/uniprot/{}.fasta".format(protein)
        n_attempt = 3
        attempt = 0
        while attempt < n_attempt:
            r = requests.get(address)
            if r.status_code == 200:
                gene = self.get_gene_from_fasta(r.text)
                return gene

            attempt += 1
            log.warning(f"attempt: {attempt}")
            log.warning(f"status_code: {r.status_code}")
            time.sleep(5)

        log.error(f"COULD NOT RETRIEVE GENE FOR PROTEIN: {protein}")
        return "N/A"

    @staticmethod
    def get_gene_from_fasta(fasta_text):
        info_line = fasta_text.split('\n')[0]

        # pattern = re.compile(r"GN=(.+)")
        pattern = re.compile(r"GN=(\S+)(\s)")

        # Does not exists in UniProt server.
        if re.search(pattern, info_line) is None:
            return "N/A"

        gene = re.search(pattern, info_line).group(1)

        return gene


class GeneIDFetcher:
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
            log.critical("No mapping data found. Make sure you have the mapping data.\n"
                         f"Path: {self.MAPPING_DATA_NAME} is not valid.")
            raise FileNotFoundError("No mapping data found!")

        self.mapping_data = mapping_data
        log.info("Mapping data loaded.")

    def fetch(self, protein, suppress_warning=True):
        log.debug(f"Fetching gene id of protein {protein}")

        try:
            manual_mapping = {
                "Q53F85": "DAXX",
                "B4DWA2": "CASP7",
                "Q13748": "TUBA3C",
            }
            if protein in manual_mapping:
                if not suppress_warning:
                    log.warning(f"Manual mapping used: {protein} : {manual_mapping[protein]} ")

                gene = manual_mapping[protein]

            else:
                gene = self.mapping_data.loc[protein]["GENE"]

        except KeyError:
            gene = "NAN"

        return gene

# pairs = [
#     ("Q17R98", "ZNF827"),
#     ("P24864", "CCNE1"),
#     ("Q9H0D6", "XRN2"),
#     ("O43149", "ZZEF1"),
#     ("Q9P253", "VPS18"),
#     ("Q68CR1", "SEL1L3"),
#     ("P35499", "SCN4A"),
#     ("Q96DU3", "SLAMF6"),
#     ("Q6ZS81", "WDFY4"),
#     ("Q99715", "COL12A1"),
#     ("P0C7W0", "PRR29"),
#     ("Q8NEK5", "ZNF548"),
#     ("P49755", "TMED10"),
#     # additional
#     ("Q14315", "FLNC"),
#     ("P55287", "CDH11"),
#     ("P78509", "RELN"),
#     ("Q8TC44", "POC1B"),
#     ("Q9NVI1", "FANCI"),
#     ("Q8WX93", "PALLD"),
#     ("Q9P2N4", "ADAMTS9"),
#     ("Q13740", "ALCAM"),
#     ("P98161", "PKD1"),
#     ("Q9Y255", "PRELID1"),
# ]
#
# gene_id_retriever = GeneIDRetriever()
#
# for pair in pairs:
#     assert gene_id_retriever.fetch(pair[0]) == pair[1], (pair[0], pair[1])
