from utils.data_loader import DataLoaderHelper
from utils.log_script import ColorHandler, MyLog

import logging

log = MyLog("debug_runner", level=logging.INFO)
log.addHandler(ColorHandler())

HIGH_CONFIDENCE_DISRUPTIVE_DATA_PATH = "../HighConfidenceDisruptiveData"


def main():
    DataLoaderHelper(HIGH_CONFIDENCE_DISRUPTIVE_DATA_PATH).scrap()


if __name__ == '__main__':
    main()
