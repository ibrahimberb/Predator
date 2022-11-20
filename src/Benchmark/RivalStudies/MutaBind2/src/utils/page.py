from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.common.exceptions import TimeoutException
from selenium.common.exceptions import NoSuchElementException

import time
from tqdm import tqdm

from .conf import TEMP_FILES_PATH
from .misc import get_PDB_file_path, wait, create_mutation_file

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

    job_id = driver.find_element_by_xpath("/html/body/div[2]/div[2]/div[1]/div/ul/li[1]").text
    pdb_id_registered = driver.find_element_by_xpath("/html/body/div[2]/div[2]/div[1]/div/ul/li[1]").text
    mutated_chain = driver.find_element_by_xpath("/html/body/div[2]/div[2]/div[2]/div/table/tbody/tr/td[1]").text
    mutation_registered = driver.find_element_by_xpath("/html/body/div[2]/div[2]/div[2]/div/table/tbody/tr/td[2]").text
    ddG_predicted = driver.find_element_by_xpath("/html/body/div[2]/div[2]/div[2]/div/table/tbody/tr/td[3]").text
    interface = driver.find_element_by_xpath("/html/body/div[2]/div[2]/div[2]/div/table/tbody/tr/td[4]").text
    deleterious = driver.find_element_by_xpath("/html/body/div[2]/div[2]/div[2]/div/table/tbody/tr/td[5]").text

    log.debug(f" > job_id: {job_id}")
    log.debug(f" > pdb_id_registered: {pdb_id_registered}")
    log.debug(f" > mutated_chain: {mutated_chain}")
    log.debug(f" > mutation_registered: {mutation_registered}")
    log.debug(f" > ddG_predicted: {ddG_predicted}")
    log.debug(f" > interface: {interface}")
    log.debug(f" > deleterious: {deleterious}")


    # chain = chain.split()[1]
    # position = position.split()[1]
    # wildtype = wildtype.split()[1]
    # mutant = mutant.split()[1]
    # dis_closest_partner = dis_closest_partner.split()[-1]

    mutation_details = {
        "job_id": job_id,
        "pdb_id_registered": pdb_id_registered,
        "mutated_chain": mutated_chain,
        "mutation_registered": mutation_registered,
        "ddG_predicted": ddG_predicted,
        "interface": interface,
        "deleterious": deleterious,
    }

    return mutation_details


def get_url(driver):
    url = driver.current_url
    log.debug("current_url: {}".format(url))
    return url


def upload_pdb_file(driver, pdb_id):
    # Locate the PDB file
    upload_file_path = get_PDB_file_path(pdb_id)
    # Upload the file.
    upload_file_button = driver.find_element_by_id("pdb_file_input")
    upload_file_button.send_keys(upload_file_path)
    log.debug('Uploaded PDB file: {}'.format(upload_file_path))


def upload_mutation_file(driver, mutation_filepath):
    upload_file_button = driver.find_element_by_name("muta_file")
    upload_file_button.send_keys(mutation_filepath)
    log.debug('Uploaded mutation file: {}'.format(mutation_filepath))



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


def click_single_mutations(driver):
    WebDriverWait(driver, 3)
    click_button_by_id(driver, "single")
    log.debug(f"Single mutation is clicked.")


def select_mutations(driver, protein, chain_to_mutate, mutation):
    """
    # It seems easier to upload a file.
    # The example file looks like this;
    # #Chain to Mutate	Residue	Mutant Residue
    # #split by tab
    # E	A18	G
    # E	G7	L
    # I	Q20	F
    # E	G7	A
    """
    # click on "Upload file" first.
    button_wait = WebDriverWait(driver, 15)
    button_upload_file_xpath = '/html/body/div[2]/div[2]/div[2]/ul/li[2]/a'
    button_wait.until(EC.visibility_of_element_located((By.XPATH, button_upload_file_xpath)))
    click_button_by_xpath(driver, button_upload_file_xpath)
    wait(1)

    # create and upload mutation txt file.
    filename = f"{protein}_{mutation}.txt"
    mutation_filepath = create_mutation_file(
        path=TEMP_FILES_PATH, filename=filename, mutation=mutation, chain_to_mutate=chain_to_mutate
    )
    wait(1)
    upload_mutation_file(driver, mutation_filepath)

    # Click on "SUBMIT JOB".
    wait(1)
    click_button_by_id(driver, "identify_submit_button")
    log.debug("Clicking on Submit Job.")


def select_partners_of_interactions(driver, chain_1, chain_2, allowed_timeout=300):
    try:
        find_chains_wait = WebDriverWait(driver, allowed_timeout)
        find_chains_wait.until(EC.visibility_of_element_located((By.ID, "gallery")))
        # chains = driver.find_elements_by_id("gallery")
        chains = driver.find_elements_by_xpath('//*[@id="gallery"]/li')
        log.info(f"Number of chains: {len(chains)}")

        for item in chains:
            log.info(f"Chain: {item.text.split()}")

        if len(chains) < 2:
            log.critical("There are less than two chains..")
            raise ValueError("Number of chains are not enough.")

        chain_to_loc = {}
        for i, item in enumerate(chains):
            chain_on_page = item.text.split()[0]
            if chain_1 == chain_on_page:
                chain_to_loc[chain_1] = i

            if chain_2 == chain_on_page:
                chain_to_loc[chain_2] = i

        if len(chain_to_loc) != 2:
            log.critical(f"Locations are not quite determined: {chain_to_loc}")
            raise ValueError("We need two chains, I guess?")

        log.info(f"CHAIN TO LOCATIONS: {chain_to_loc}")
        # Select the first chain
        click_button_by_xpath(
            driver, f'/html/body/div[2]/div[2]/div[4]/div[1]/div/div/ul/li[{chain_to_loc[chain_1] + 1}]/a[1]'
        )

        # Select the second chain
        click_button_by_xpath(
            driver, f'/html/body/div[2]/div[2]/div[4]/div[1]/div/div/ul/li[{chain_to_loc[chain_2] + 1}]/a[2]'
        )

        # Click on `next` button.
        wait(10)
        log.debug("Clicking on `NEXT`")
        click_button_by_xpath(driver, '/html/body/div[2]/div[2]/div[4]/div[2]/button')

    except TimeoutException:
        raise


def click_button_by_id(driver, element_id, allowed_timeout=5):
    try:
        click_button_wait = WebDriverWait(driver, allowed_timeout)
        click_button_wait.until(EC.visibility_of_element_located((By.ID, element_id)))
        driver.find_element_by_id(element_id).click()
    except TimeoutException:
        raise


def click_button_by_xpath(driver, element_xpath, allowed_timeout=25):
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
    process_submission_allowed_time = 2 * 60 * 60  # 2 hrs
    wait_results_page_load = WebDriverWait(driver, 20)  # 20 seconds to load
    submission_text_xpath = '/html/body/div[2]/div[2]/h2'
    try:
        wait_results_page_load.until(EC.visibility_of_element_located((By.XPATH, submission_text_xpath)))
        results_xpath = '/html/body/div[2]/div[2]/h2'
        if driver.find_element_by_xpath(results_xpath).text in ["Job submitted..", "Calculating.."]:
            log.debug(driver.find_element_by_xpath(results_xpath).text)
            computation_time_elapsed_seconds = 0
            job_completed = False
            check_cooldown = 20  # check in every n seconds.

            pbar = tqdm(total=process_submission_allowed_time, desc='still not completed. waiting', position=0,
                        leave=True)
            while computation_time_elapsed_seconds < process_submission_allowed_time and not job_completed:
                computation_time_elapsed_seconds += check_cooldown
                time.sleep(check_cooldown)

                wait_results_page_load.until(EC.visibility_of_element_located((By.XPATH, results_xpath)))
                if driver.find_element_by_xpath(results_xpath).text.strip() == "Results":
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
    WebDriverWait(driver, 10)

    # Finding message "Error we have the following requirements for the input file,
    # and your submitted file does not meet one or more."

    validating_text_xpath = "/html/body/div[2]/div[2]/div[1]/div[1]"
    if check_exists_by_xpath(driver, validating_text_xpath):
        message = driver.find_element_by_xpath(validating_text_xpath).text
        if message.split()[0].lower() == "error:":
            wait(1)
            log.warning(message)
            return message

        job_submitted_text_xpath = "/html/body/div[2]/div[2]/h2"
        message = driver.find_element_by_xpath(job_submitted_text_xpath).text
        wait(20)
        if check_exists_by_xpath(driver, job_submitted_text_xpath) and message in ["Job submitted..", "Calculating.."]:
            log.debug("Validation success! We are processing your submission ...")

        else:
            raise Exception("Unknown error occurred.")

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
