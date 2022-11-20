from selenium import webdriver
from .conf import TEMP_DOWNLOAD_FOLDER_PATH, DRIVER_PATH, HEADLESS, PATH_TO_DEV_NULL


def initialize_driver():
    # Set the Options.
    profile = webdriver.FirefoxProfile()
    profile.set_preference('browser.download.folderList', 2)  # custom location
    profile.set_preference('browser.download.manager.showWhenStarting', False)
    profile.set_preference("browser.download.dir", str(TEMP_DOWNLOAD_FOLDER_PATH))
    profile.set_preference('browser.helperApps.neverAsk.saveToDisk', 'text/plain')  # type of file to download

    options = webdriver.FirefoxOptions()
    options.headless = HEADLESS

    # Set the driver
    driver = webdriver.Firefox(
        options=options,
        executable_path=DRIVER_PATH,
        firefox_profile=profile,
        service_log_path=PATH_TO_DEV_NULL
    )
    driver.set_window_position(1024, 0)
    driver.set_window_size(1920 - 1024, 1040)

    return driver


def close_driver(driver):
    driver.quit()
