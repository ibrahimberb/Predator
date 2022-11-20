from src.utils.errors import MissionAbortError
from src.utils.record import Record, Entry, post_check_entry_data
from utils.driver_conf import initialize_driver, close_driver
from utils.log_script import ColorHandler
from utils.resources import terminate_firefox_processes
from utils.conf import MUTABIND2_URL
from utils.page import (
    upload_pdb_file,
    fill_mutation,
    fill_pdb_accession,
    fill_chain,
    submit_input,
    process_validate_input_data,
    process_submission,
    SubmissionResponseMessages,
    retrieve_completed_results,
    get_url, get_mutation_details, click_single_mutations, select_partners_of_interactions, select_mutations
)

import logging

log = logging.Logger("debug_runner", level=logging.DEBUG)
log.addHandler(ColorHandler())


class MyMutaBind2API:
    def __init__(
            self,
            mutation_effect_label,
            protein,
            mutation,
            interactor,
            pdb_id,
            chain_id_1,
            chain_id_2,
            upload_our_pdb
    ):
        log.info("\nMutaBind2 API is started")
        self.mutation_effect_label = mutation_effect_label
        self.protein = protein
        self.mutation = mutation
        self.interactor = interactor
        self.pdb_id = pdb_id
        self.chain_id_1 = chain_id_1
        self.chain_id_2 = chain_id_2
        self.upload_our_pdb = upload_our_pdb
        self.mutation_details = {}

        log.debug(f"MUTATION EFFECT: {self.mutation_effect_label}")
        log.debug(f"PROTEIN: {self.protein}")
        log.debug(f"MUTATION: {self.mutation}")
        log.debug(f"INTERACTOR: {self.interactor}")
        log.debug(f"PDB_ID: {self.pdb_id}")
        log.debug(f"CHAIN_ID_1: {self.chain_id_1}")
        log.debug(f"CHAIN_ID_2: {self.chain_id_2}")

        self.driver = None
        self.url = None

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
        self.driver.get(MUTABIND2_URL)

    def save_failed_results(self, error):
        entry = Entry(
            mutation_effect_label=self.mutation_effect_label,
            protein=self.protein,
            mutation=self.mutation,
            interactor=self.interactor,
            pdb_id=self.pdb_id,
            chain_id_1=self.chain_id_1,
            chain_id_2=self.chain_id_2,
            submitted=1,
            saved=0,
            error_encountered=error,
            result_url=None,
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
            chain_id_1=self.chain_id_1,
            chain_id_2=self.chain_id_2,
            submitted=1,
            saved=1,
            error_encountered=None,
            result_url=self.url,
            mutation_details=self.mutation_details
        )

        self.record.add_entry(entry, successful_entry=True)

    def log_results(self):
        print(f"{self.protein=}")
        print(f"{self.mutation=}")
        print(f"{self.interactor=}")
        print(f"{self.pdb_id=}")
        print(f"{self.chain_id_1=}")
        print(f"{self.chain_id_2=}")
        print(f"{self.url=}")

    def handle_input(self):
        # Step 1 - Select Protein Complex
        if self.upload_our_pdb:
            upload_pdb_file(self.driver, self.pdb_id)
        else:
            fill_pdb_accession(self.driver, self.pdb_id)

        # click on single mutation
        click_single_mutations(self.driver)

        # Step 2 - Select Protein Complex
        select_partners_of_interactions(self.driver, self.chain_id_1, self.chain_id_2)

        # Step 3 - Select Mutations
        select_mutations(
            driver=self.driver,
            protein=self.protein,
            chain_to_mutate=self.chain_id_1,
            mutation=self.mutation
        )

        error_message = process_validate_input_data(self.driver)
        if error_message is not None:
            self.save_failed_results(error_message)
            close_driver(self.driver)
            raise MissionAbortError

        self.url = get_url(self.driver)

    def handle_process_completed(self):
        log.info("Completed successfully.")
        self.mutation_details = get_mutation_details(self.driver)
