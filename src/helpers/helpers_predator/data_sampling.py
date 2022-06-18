import logging
import pandas as pd
import numpy as np

from .common import print_annotation

logger = logging.getLogger(__name__)


def prepare_data_spsm(data: pd.DataFrame, random_seed):
    """
    Prepares data with Single Protein Single Mutation form.
    :param random_seed:
    :param data:
    :return: sampled_train_data
    """

    # Control the behavior of randomization in a reproducable manner.
    np.random.seed(random_seed)

    # Get the unique proteins from `UniProt_ID` column.
    unique_proteins = list(data['UniProt_ID'].unique())  # todo: maybe wrap with `sorted`.

    logger.debug('Number of `unique_proteins`:', len(unique_proteins))
    logger.debug('First five proteins: {}'.format(unique_proteins[:5]))

    sampled_row_dataframes = []
    for unique_protein in unique_proteins:
        sampled_row_dataframes.append(data[data['UniProt_ID'] == unique_protein].sample())

    # Merge row dataframes into single dataframe, stack rows on top of each other.
    sampled_train_data = pd.concat(sampled_row_dataframes, axis='rows')

    # Reset index of the dataframe to avoid any possible errors
    sampled_train_data.reset_index(drop=True, inplace=True)

    logger.debug(f"Dimensions of sampled_dataframe: {sampled_train_data.shape}")

    return sampled_train_data


# FIXME
def prepare_data_spmm(data: pd.DataFrame):
    """
    Prepares data with Single Protein Multiple Mutation form.
    TODO: docstring
    :return: sampled_train_data
    """

    # Introducing new column `Protein_Mutation`, containing (protein, mutation) tuple.
    data['Protein_Mutation'] = data.apply(lambda x: (x['UniProt_ID'], x['Mutation']), axis=1)

    # Get the unique (protein, mutation) pairs from `Protein_Mutation` column.
    unique_protein_mutation_pairs = list(data['Protein_Mutation'].unique())

    # Number of unique (protein, mutation) pairs.
    print('Number of `unique_protein_mutation_pairs`:', len(unique_protein_mutation_pairs))

    # First five (protein, mutation) pairs
    print(unique_protein_mutation_pairs[:5])

    sampled_row_dataframes = []
    for unique_protein_mutation in unique_protein_mutation_pairs:
        sampled_row_dataframes.append(
            data[data['Protein_Mutation'] == unique_protein_mutation].sample())

    # Merge row dataframes into single dataframe, stack rows on top of each other.
    sampled_train_data = pd.concat(sampled_row_dataframes)

    # Reset index of the dataframe to avoid any possible errors
    sampled_train_data.reset_index(drop=True, inplace=True)

    # Dimensions of dataframe
    print_annotation(f"Dimensions of sampled_dataframe: {sampled_train_data.shape}")

    return sampled_train_data
