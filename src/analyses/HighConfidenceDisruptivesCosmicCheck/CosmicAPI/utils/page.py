import re
import time

from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.common.exceptions import TimeoutException
from selenium.common.exceptions import NoSuchElementException

import pandas as pd

from .conf import USER_CREDENTIALS_PATH

from .misc import load_credentials

from .log_script import ColorHandler, MyLog
import logging

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# logging.basicConfig(level=logging.INFO, format='%(message)s')

log = MyLog("debug_runner", level=logging.INFO)
log.addHandler(ColorHandler())


class CosmicTierColors:
    TIER_COLOR_MAPPING = {
        ('143', '25', '27'): "TIER_1",
        ('211', '87', '39'): "TIER_2",
        ('218', '217', '130'): "TIER_3",
        ('253', '251', '250'): "OTHER_MUTATION",
    }

    def __init__(self):
        pass

    def get_tier(self, style):
        rgb_tuple = re.search(r"rgb\((\d+), (\d+), (\d+)\)", style).groups()
        return self.TIER_COLOR_MAPPING[rgb_tuple]


def enter_credentials(driver):
    log.debug("Filling the credentials..")
    credentials = load_credentials(USER_CREDENTIALS_PATH)
    input_email = driver.find_element_by_id("id_username")
    input_email.send_keys(credentials["email"])
    log.debug("Email entered.")

    input_password = driver.find_element_by_id("id_password")
    input_password.send_keys(credentials["password"])
    log.debug("Password entered.")


def login(driver):
    log.debug("Logging in..")
    enter_credentials(driver)
    submit_button_xpath = '/html/body/div/div/form/input[2]'
    click_button_by_xpath(driver, submit_button_xpath)
    log.debug("Logged in successfully.")


def retrieve(driver):
    no_results_found_xpath = '//*[@id="root"]/div[2]/div/div[2]/div/h1'
    codon_intro_xpath = '//*[@id="root"]/div[2]/div/div[2]/div[1]'

    max_attempt = 10
    n_attempt = 0
    while True:

        try:
            if check_exists_by_xpath(driver, no_results_found_xpath):
                log.critical("No results found!")
                error_message_xpath = '//*[@id="root"]/div[2]/div/div[2]/div/div/h2'
                error_message = driver.find_element_by_xpath(error_message_xpath).text
                log.error(error_message)
                error_info = {"error_message": error_message}
                return (
                    SubmissionResponseMessages.NoResultsFound,
                    error_info
                )

            elif check_exists_by_xpath(driver, codon_intro_xpath):
                log.success("Gene is found in COSMIC.")
                grid_overlay_info = get_grid_overlay_info(driver)
                return (
                    SubmissionResponseMessages.ResultsFound,
                    grid_overlay_info
                )

        except NoSuchElementException:
            if n_attempt < max_attempt:
                log.error(f"Something's not right, I'll try again .. (n_attempt = {n_attempt})")
                time.sleep(10)
                n_attempt += 1

            else:
                log.critical(f"Failed to load necessary items.")
                raise




def get_grid_overlay_info(driver):
    grid_overlays = driver.find_elements_by_class_name("grid-overlay")
    codon_xpath = '/html/body/div/div[2]/div/div[1]/div[1]/div[2]/span[2]/span'

    try:
        gene_xpath = '/html/body/div/div[2]/div/div[1]/div[1]/div[2]/span[1]/a'
        gene = driver.find_element_by_xpath(gene_xpath).text
    except NoSuchElementException:
        gene_xpath = '/html/body/div/div[2]/div/div[1]/div[1]/div[2]/span[1]/p'
        gene = driver.find_element_by_xpath(gene_xpath).text

    codon = driver.find_element_by_xpath(codon_xpath).text

    try:
        CGC_status_xpath = '/html/body/div/div[2]/div/div[1]/div[1]/div[2]/span[3]/a'
        CGC_status = driver.find_element_by_xpath(CGC_status_xpath).text
    except NoSuchElementException:
        CGC_status_xpath = '//*[@id="root"]/div[2]/div/div[1]/div[1]/div[2]/span[3]/p'
        CGC_status = driver.find_element_by_xpath(CGC_status_xpath).text

    log.debug(f"number of grid_overlays: {len(grid_overlays)}")
    log.debug("- - - - - - - - - ")
    grid_overlay_entries = []
    for grid_overlay in grid_overlays:
        # log.debug(f"{grid_overlay.text}")
        log.debug(CosmicTierColors().get_tier(grid_overlay.get_attribute("style")))
        detailed_info = get_detailed_info(
            CosmicTierColors().get_tier(grid_overlay.get_attribute("style")),
            grid_overlay.text,
        )
        grid_overlay_entries.append(detailed_info)
        log.debug(f"detailed_info: {detailed_info}")
        log.debug("- - - - - - - - - ")

    grid_overlay_data = pd.DataFrame(
        grid_overlay_entries,
        columns=[
            "TIER", "COSMIC_MUT", "NUCLEOTIDE_CHANGE", "AA_CHANGE", "MUTATION_TYPE"
        ]
    )

    grid_overlay_data.insert(loc=0, column="GENE", value=gene)
    grid_overlay_data.insert(loc=1, column="CODON", value=codon)
    grid_overlay_data.insert(loc=2, column="CGC_STATUS", value=CGC_status)

    # Assuming codons (positions) ranked by TIER, the first one should be the most significant tier.
    most_significant_codon_tier = grid_overlay_data["TIER"][0]

    output = {
        "grid_overlay_data": grid_overlay_data,
        "gene": gene,
        "codon": codon,
        "CGC_status": CGC_status,
        "most_significant_codon_tier": most_significant_codon_tier,
    }

    return output


def get_detailed_info(*args):
    info = []
    for item in args:
        info.extend(item.split("\n"))

    return tuple(info)


class SubmissionResponseMessages:
    NoResultsFound = "NoResultsFound"
    ResultsFound = "ResultsFound"


def click_button_by_xpath(driver, element_xpath, allowed_timeout=5):
    try:
        click_button_wait = WebDriverWait(driver, allowed_timeout)
        click_button_wait.until(EC.visibility_of_element_located((By.XPATH, element_xpath)))
        driver.find_element_by_xpath(element_xpath).click()
    except TimeoutException:
        raise


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
