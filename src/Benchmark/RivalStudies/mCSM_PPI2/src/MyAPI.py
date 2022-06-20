from src.utils.errors import MissionAbortError
from src.utils.record import Record, Entry, post_check_entry_data
from utils.driver_conf import initialize_driver, close_driver
from utils.log_script import ColorHandler
from utils.resources import terminate_firefox_processes
from utils.conf import MCSM_PPI_MANY_URL
from utils.page import (
    upload_file,
    fill_mutation,
    fill_pdb_accession,
    fill_chain,
    submit_input,
    process_validate_input_data,
    process_submission,
    SubmissionResponseMessages,
    retrieve_completed_results,
    get_url, get_mutation_details
)

import logging

log = logging.Logger("debug_runner", level=logging.DEBUG)
log.addHandler(ColorHandler())


class MyMCSMPP2API:
    def __init__(
            self,
            mutation_effect_label,
            protein,
            mutation,
            interactor,
            pdb_id,
            chain_id,
            upload_our_pdb
    ):
        log.info("\nmCSM-PP2 API is started")
        self.mutation_effect_label = mutation_effect_label
        self.protein = protein
        self.mutation = mutation
        self.interactor = interactor
        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.upload_our_pdb = upload_our_pdb
        self.mutation_details = {}

        self.driver = None
        self.url = None
        self.predicted_value = None
        self.predicted_change = None

        self.record = Record()

        self.open_default_page()
        self.handle_input()

        response = process_submission(self.driver)
        if response == SubmissionResponseMessages.COMPLETED:
            self.handle_process_completed()
            self.save_successful_results()
            close_driver(self.driver)
        else:
            raise

        self.log_results()
        # terminate_firefox_processes()
        # # # the end # # #

    def open_default_page(self):
        self.driver = initialize_driver()
        self.driver.get(MCSM_PPI_MANY_URL)

    def save_failed_results(self, error):
        entry = Entry(
            mutation_effect_label=self.mutation_effect_label,
            protein=self.protein,
            mutation=self.mutation,
            interactor=self.interactor,
            pdb_id=self.pdb_id,
            chain_id=self.chain_id,
            submitted=1,
            saved=0,
            error_encountered=error,
            result_url=None,
            predicted_affinity_change_value=None,
            predicted_affinity_change=None,
            mutation_details=None
        )

        self.record.add_entry(entry, successful_entry=False)

    def save_successful_results(self):
        entry = Entry(
            mutation_effect_label=self.mutation_effect_label,
            protein=self.protein,
            mutation=self.mutation,
            interactor=self.interactor,
            pdb_id=self.pdb_id,
            chain_id=self.chain_id,
            submitted=1,
            saved=1,
            error_encountered=None,
            result_url=self.url,
            predicted_affinity_change_value=self.predicted_value,
            predicted_affinity_change=self.predicted_change,
            mutation_details=self.mutation_details
        )

        self.record.add_entry(entry, successful_entry=True)

    def log_results(self):
        print(f"{self.protein=}")
        print(f"{self.mutation=}")
        print(f"{self.interactor=}")
        print(f"{self.pdb_id=}")
        print(f"{self.chain_id=}")
        print(f"{self.url=}")
        print(f"{self.predicted_value=}")
        print(f"{self.predicted_change=}")

    def handle_input(self):
        if self.upload_our_pdb:
            upload_file(self.driver, self.pdb_id)
        else:
            fill_pdb_accession(self.driver, self.pdb_id)

        fill_mutation(self.driver, self.mutation)
        fill_chain(self.driver, self.chain_id)
        submit_input(self.driver)
        error_message = process_validate_input_data(self.driver)
        if error_message is not None:
            self.save_failed_results(error_message)
            close_driver(self.driver)
            raise MissionAbortError

        self.url = get_url(self.driver)

    def handle_process_completed(self):
        log.info("Completed successfully.")
        predicted_value, predicted_change = retrieve_completed_results(self.driver)
        self.predicted_value = predicted_value
        self.predicted_change = predicted_change
        self.mutation_details = get_mutation_details(self.driver)
