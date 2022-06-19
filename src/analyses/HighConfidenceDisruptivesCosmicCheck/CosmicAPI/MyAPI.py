from utils.driver_conf import initialize_driver

import time

from utils.log_script import ColorHandler, MyLog
from utils.misc import get_cosmic_url, save_CGC_data, get_residue_position
from utils.page import login, retrieve, SubmissionResponseMessages

import logging

from utils.record import Record, Entry

log = MyLog("debug_runner", level=logging.INFO)
log.addHandler(ColorHandler())


class MyCosmicAPI:
    def __init__(
            self,
            gene,
            mutation,
    ):
        log.info("=========================================")
        log.info("Accessing COSMIC Cancer Mutation Census ..")
        log.info(f"GENE: {gene}")
        log.info(f"MUTATION: {mutation}")
        self.gene = gene
        self.residue_position = get_residue_position(mutation)
        self.url = get_cosmic_url(gene, mutation)
        self.record = Record()
        self.CGC_status = None
        self.most_significant_codon_tier = None

        self.driver = None
        self.resulting_data = None

        self.run()

    def query(self):
        log.info("Querying ..")
        self.driver = initialize_driver()
        self.driver.get(self.url)
        time.sleep(3)
        login(self.driver)
        time.sleep(3)

        output = retrieve(self.driver)
        return output

    def run(self):
        # if already recorded, then no need to open browser.
        if self.record.is_already_recorded(
                gene=self.gene, residue_position=self.residue_position
        ):
            log.info("Already recorded.")
            return

        # Accessing through Web API
        response, output = self.query()

        if response == SubmissionResponseMessages.ResultsFound:
            self.resulting_data = output["grid_overlay_data"]
            self.CGC_status = output["CGC_status"]
            self.most_significant_codon_tier = output["most_significant_codon_tier"]
            # Saving CGC detailed information.
            save_CGC_data(self.resulting_data, self.gene, self.residue_position)
            submitted = 1
            downloaded = 1
            error_encountered = None

        elif response == SubmissionResponseMessages.NoResultsFound:
            submitted = 1
            downloaded = 0
            error_encountered = output["error_message"]

        else:
            raise

        log.red(f"CGC_status: {self.CGC_status}")

        # Saving the record.
        entry = Entry(
            gene=self.gene,
            residue_position=self.residue_position,
            submitted=submitted,
            downloaded=downloaded,
            error_encountered=error_encountered,
            most_significant_codon_tier=self.most_significant_codon_tier,
            CGC_status=self.CGC_status,
            url=self.url
        )
        self.record.add_entry(entry)
        self.driver.quit()
        return


# # URL = r"https://cancer.sanger.ac.uk/cmc/gene/TP53/codon/389"
# input_gene = "TP53"
# input_gene_2 = "TP53"
# input_gene_3 = "TP53"
# input_gene_4 = "NOSUCHGENE"
# # input_mutation = "X389X"
# # input_mutation = "X304X"
# input_mutation = "X220X"
# input_mutation_2 = "X304X"
# input_mutation_3 = "X123456789X"
# input_mutation_4 = "X304X"

# https://cancer.sanger.ac.uk/cmc/gene/TP53/codon/220

# CosmicScraper(gene=input_gene, mutation=input_mutation)
# CosmicScraper(gene=input_gene_2, mutation=input_mutation_2)
# CosmicScraper(gene=input_gene_3, mutation=input_mutation_3)
# CosmicScraper(gene=input_gene_4, mutation=input_mutation_4)
