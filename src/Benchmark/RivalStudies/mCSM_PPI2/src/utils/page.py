from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.common.exceptions import TimeoutException
from selenium.common.exceptions import NoSuchElementException

import time
from tqdm import tqdm

from .errors import ValidationProcessFailedError
from .misc import get_PDB_file_path, wait

from .log_script import ColorHandler
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')

log = logging.Logger("debug_runner", level=logging.DEBUG)
log.addHandler(ColorHandler())


class SubmissionResponseMessages:
    STILL_PROCESSING = "STILL_PROCESSING"
    COMPLETED = "COMPLETED"
    RESULT_PAGE_NOT_LOADED = "RESULT_PAGE_NOT_LOADED"


def get_mutation_details(driver):
    chain = driver.find_element_by_xpath("/html/body/main/div[2]/div/div[1]/div/div/div/div/div[2]/p[2]").text
    position = driver.find_element_by_xpath("/html/body/main/div[2]/div/div[1]/div/div/div/div/div[2]/p[3]").text
    wildtype = driver.find_element_by_xpath("/html/body/main/div[2]/div/div[1]/div/div/div/div/div[2]/p[4]").text
    mutant = driver.find_element_by_xpath("/html/body/main/div[2]/div/div[1]/div/div/div/div/div[2]/p[5]").text
    dis_closest_partner = driver.find_element_by_xpath("/html/body/main/div[2]/div/div[1]/div/div/div/div/div[2]/p[6]").text

    chain = chain.split()[1]
    position = position.split()[1]
    wildtype = wildtype.split()[1]
    mutant = mutant.split()[1]
    dis_closest_partner = dis_closest_partner.split()[-1]

    log.debug(f" > chain: {chain}")
    log.debug(f" > position: {position}")
    log.debug(f" > wildtype: {wildtype}")
    log.debug(f" > mutant: {mutant}")
    log.debug(f" > dis_closest_partner: {dis_closest_partner}")

    mutation_details = {
        "chain": chain,
        "position": position,
        "wildtype": wildtype,
        "mutant": mutant,
        "dis_closest_partner": dis_closest_partner,
    }

    return mutation_details


def get_url(driver):
    url = driver.current_url
    log.debug("current_url: {}".format(url))
    return url


def upload_file(driver, pdb_id):
    # Locate the PDB file
    upload_file_path = get_PDB_file_path(pdb_id)
    # Upload the file.
    upload_file_button = driver.find_element_by_id("pdb_file_single")
    upload_file_button.send_keys(upload_file_path)
    log.debug('Uploaded file: {}'.format(upload_file_path))


def fill_pdb_accession(driver, pdb_association):
    input_element_pdb_accession = driver.find_element_by_id("pdb_accession_single")
    input_element_pdb_accession.send_keys(pdb_association.upper())
    log.debug(f"PDB association {pdb_association} is entered.")


def fill_mutation(driver, mutation):
    mutation = mutation.upper()
    input_element_mutation_single = driver.find_element_by_id("mutation_single")
    input_element_mutation_single.send_keys(mutation)
    log.debug(f"Mutation {mutation} is entered.")


def fill_chain(driver, chain):
    input_element_chain_single = driver.find_element_by_id("chain_single")
    input_element_chain_single.send_keys(chain)
    log.debug(f"Chain {chain} is entered.")


def submit_input(driver):
    WebDriverWait(driver, 3)
    submit_button_xpath = '//*[@id="singlePredictionForm"]/div/div[2]/div/div[1]/button'
    click_button_by_xpath(driver, submit_button_xpath)
    log.debug(f"Submitting the input ..")


def click_button_by_xpath(driver, element_xpath, allowed_timeout=5):
    try:
        click_button_wait = WebDriverWait(driver, allowed_timeout)
        click_button_wait.until(EC.visibility_of_element_located((By.XPATH, element_xpath)))
        driver.find_element_by_xpath(element_xpath).click()
    except TimeoutException:
        raise


def retrieve_completed_results(driver):
    predicted_value_xpath = '//*[@id="ppi2Prediction"]'
    predicted_change_xpath = "/html/body/main/div[2]/div/div[1]/div/div/div/div/div[1]/p[2]/i"

    predicted_value = driver.find_element_by_xpath(predicted_value_xpath).text
    predicted_value = float(predicted_value.split(" ")[0])
    predicted_change = driver.find_element_by_xpath(predicted_change_xpath).text[1:-1]

    log.debug("Results:")
    log.debug(f"predicted_value: {predicted_value}")
    log.debug(f"predicted_change: {predicted_change}")

    return predicted_value, predicted_change


def process_submission(driver):
    log.info("Computing input")
    process_submission_allowed_time = 60 * 10 * 3
    wait_results_page_load = WebDriverWait(driver, 10)  # 10 seconds to load
    submission_text_xpath = '//*[@id="index-banner"]/div/div/div/h4'
    try:
        wait_results_page_load.until(EC.visibility_of_element_located((By.XPATH, submission_text_xpath)))
        results_xpath = '/html/body/main/div[2]/div/div[1]/div/div/div/span'
        if not check_exists_by_xpath(driver, results_xpath):
            computation_time_elasped_seconds = 0
            job_completed = False
            check_cooldown = 1  # check in every 1 seconds.

            pbar = tqdm(total=process_submission_allowed_time, desc='still not completed. waiting', position=0,
                        leave=True)
            while computation_time_elasped_seconds < process_submission_allowed_time and not job_completed:
                computation_time_elasped_seconds += check_cooldown
                time.sleep(check_cooldown)
                if check_exists_by_xpath(driver, results_xpath):
                    job_completed = True

                pbar.update(check_cooldown)
            pbar.close()

            if not job_completed:
                # logging.info("I had enough.. can't wait any longer.")
                return SubmissionResponseMessages.STILL_PROCESSING

            else:
                # job completed
                return SubmissionResponseMessages.COMPLETED

        else:
            raise Exception('Unexpected error.')

    except TimeoutException:
        log.warning("Result page not loaded!")
        return SubmissionResponseMessages.RESULT_PAGE_NOT_LOADED


def process_validate_input_data(driver):
    log.debug("validating input ..")
    wait_process_page_load = WebDriverWait(driver, 1)

    # Finding message "We are checking your input"
    validating_text_xpath = "/html/body/div[3]/div[2]/h2"
    # wait_process_page_load.until(EC.visibility_of_element_located((By.XPATH, validating_text_xpath)))

    # log.debug(f"we found it {driver.find_element_by_xpath(validating_text_xpath).text}")

    if check_exists_by_xpath(driver, validating_text_xpath):
        computation_time_elasped_seconds = 0
        validation_completed = False
        check_cooldown = 1  # check in every 1 seconds.

        while computation_time_elasped_seconds < 120 and not validation_completed:
            computation_time_elasped_seconds += check_cooldown
            time.sleep(check_cooldown)
            if not check_exists_by_xpath(driver, validating_text_xpath):  # If we don't see it, it is done.
                validation_completed = True

        if not validation_completed:
            raise ValidationProcessFailedError

    processing_input_xpath = '/html/body/main/div[2]/div/div/div/h5'
    submission_error_xpath = '/html/body/main/div[1]/div/div/div/h4'

    wait(0.9)

    if check_exists_by_xpath(driver, processing_input_xpath):
        log.debug("Validation success! We are processing your submission ...")

    elif check_exists_by_xpath(driver, submission_error_xpath):
        log.critical("Submission error!")
        error_message_1 = driver.find_element_by_xpath("/html/body/main/div[2]/div/div/div[1]/p[1]").text
        error_message_2 = driver.find_element_by_xpath("/html/body/main/div[2]/div/div/div[1]/p[2]").text
        error_message_3 = driver.find_element_by_xpath("/html/body/main/div[2]/div/div/div[1]/p[3]").text
        error_message = "\n".join(
            [error_message_1, error_message_2, error_message_3]
        )
        log.warning(error_message)
        return error_message

    else:
        raise Exception("Unknown error occurred.")


def check_exists_by_id(driver, element_id):
    try:
        driver.find_element_by_id(element_id)
    except NoSuchElementException:
        return False
    return True


def check_exists_by_xpath(driver, element_xpath):
    try:
        driver.find_element_by_xpath(element_xpath)
    except NoSuchElementException:
        return False
    return True
