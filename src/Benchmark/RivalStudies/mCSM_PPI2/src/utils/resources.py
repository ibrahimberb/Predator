from .conf import ALLOWED_RAM_PERCENTAGE
from .log_script import ColorHandler

import psutil
import logging
import subprocess

log = logging.Logger('computation_utils', level=logging.WARNING)
log.addHandler(ColorHandler())


def get_ram_usage_percentage():
    return psutil.virtual_memory().percent
    # percentage of available memory
    # return psutil.virtual_memory().available * 100 / psutil.virtual_memory().total


def get_firefox_processes():
    result = subprocess.run("tasklist", shell=True, check=True, stdout=subprocess.PIPE)

    firefox_processes = [s for s in result.stdout.decode('utf-8').split()
                         if 'firefox.exe' in s]
    return firefox_processes


def terminate_firefox_processes():
    firefox_processes = get_firefox_processes()

    if firefox_processes:
        log.warning("Number of Firefox processes: {}".format(len(firefox_processes)))
        log.warning('Terminating firefox processes ..')
        subprocess.run("taskkill /F /IM firefox.exe", shell=True, check=True, stdout=subprocess.PIPE)
        log.warning('Firefox processes are terminated.')
        assert len(get_firefox_processes()) == 0

    else:
        log.warning('No Firefox process found.')
        return


def utilization_exceed() -> bool:
    # todo: chain (num_firefox_instances > 100) with OR.
    return get_ram_usage_percentage() > ALLOWED_RAM_PERCENTAGE


def check_utilization():
    if utilization_exceed():
        log.critical('{}% of RAM is occupied.'.format(get_ram_usage_percentage()))
        log.critical('RAM UTILIZATION EXCEED TO CRITICAL LEVEL!')
        raise MemoryError

    log.debug('MEMORY UTILIZATION is OK.')
    return


def handle_memory_utilization():
    if utilization_exceed():
        terminate_firefox_processes()

    check_utilization()


handle_memory_utilization()
