import pandas as pd
from .common import get_file_path
from .log_script import ColorHandler
import logging

log = logging.Logger('debug_runner', level=logging.INFO)
log.addHandler(ColorHandler())


# TODO this is class method and give a default value to project common filedir
def load_train_data(project_common_file_dir, mutations_path):
    """
    Loads the mutations dataframe, which is training data.
    :param project_common_file_dir:
    :param mutations_path:
    :return:
    """
    mutations = pd.read_csv(get_file_path(project_common_file_dir, mutations_path), sep='\t')
    log.debug(f"Size of dataframe: {mutations.shape}")
    log.debug('Dataframe head: {}'.format(mutations.head()))
    return mutations


def load_tcga_data(project_common_file_dir, tcga_path):
    """
    Loads the TCGA cancer dataset.
    :param tcga_path:
    :return:
    """
    tcga_data = pd.read_csv(get_file_path(project_common_file_dir, tcga_path), sep='\t')
    log.debug(f"Size of dataframe: {tcga_data.shape}")
    log.debug('Dataframe head: {}'.format(tcga_data.head()))
    return tcga_data
