from datetime import datetime

import pandas as pd
from scipy.stats.stats import pearsonr
from itertools import combinations
from tqdm.notebook import tqdm
from IPython.display import display
import os.path as op

from ..mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)  # INFO


def get_entries_with_protein(protein, data):
    query_data = data[
        (data["UniProt_ID"] == protein)
    ].copy()

    return query_data


def get_entries_with_pair(protein, mutation, data):
    query_data = data[
        (data["UniProt_ID"] == protein) &
        (data["Mutation"] == mutation)
        ].copy()

    return query_data


def convert_isoform(protein):
    return protein.split('-')[0]


def is_same_interactor_pair(data) -> bool:
    """
    Given data, check if interactors are same/similar.
    """
    assert data["UniProt_ID"].nunique() == 1
    assert data["Mutation"].nunique() == 1

    unique_interactors = data["Interactor_UniProt_ID"].unique()
    unique_interactors = set(map(convert_isoform, unique_interactors))

    return len(unique_interactors) == 1


def is_same_interactor(data) -> bool:
    """
    Given data, check if interactors are same/similar.
    """
    assert data["UniProt_ID"].nunique() == 1

    unique_interactors = data["Interactor_UniProt_ID"].unique()
    unique_interactors = set(map(convert_isoform, unique_interactors))

    return len(unique_interactors) == 1


def export_data_helper(
        data,
        file_name: str,
        file_extension='csv',
        overwrite=False
) -> None:
    """
    A helper function to export given data with specified name and extension.
    """

    log.debug(f"Exporting data {file_name} ..")
    file_name = op.join(file_name)
    file_date = datetime.today().strftime('%Y-%m-%d')
    file_name = f'{file_name}_{file_date}.{file_extension}'

    # Ensure the file is not exists before creating to prevent overwriting.
    if op.isfile(file_name) and not overwrite:
        log.warning(f"File {file_name} is already exist.\n"
                    "To overwrite existing file, use `overwrite=True`.")

    else:
        # Export
        data.to_csv(file_name, index=False)
        log.info(f'{file_name} is exported successfully.')


def get_corr_values(query):
    corr_score_values = []
    data_indices = query.index
    for a, b in combinations(data_indices, 2):
        # Skip the pair entries if interactor is similar/same
        if is_same_interactor(query.loc[[a, b], :]):
            log.debug(f'skip {a, b}')
            # display(query.loc[[a, b], :])
            continue

        pearson_corr_score = pearsonr(query.loc[a, :][4:], query.loc[b, :][4:])[0]
        corr_score_values.append(pearson_corr_score)

    log.debug(f"{len(corr_score_values)=}")

    try:
        corr_score = round((sum(corr_score_values) / len(corr_score_values)), 2)
    except ZeroDivisionError:
        log.critical(f"{corr_score_values=}")
        display(query)
        corr_score = "HAS ONLY ONE INTERACTOR ACROSS ALL MUTATIONS"
        # raise

    return corr_score


class MutationJustifier:
    def __init__(self, training_data_path):
        self.training_data = pd.read_csv(training_data_path)
        self.unique_proteins = sorted(
            set(self.training_data["UniProt_ID"])
        )
        self.unique_proteins_corr_data = self.get_unique_proteins_corr_data()

    def get_unique_proteins_corr_data(self):
        corr_scores = self.get_corr_scores()
        unique_proteins_corr_data = pd.DataFrame(
            self.unique_proteins, columns=["PROTEIN"]
        )
        unique_proteins_corr_data["PEARSON_CORR"] = corr_scores
        return unique_proteins_corr_data

    def get_corr_scores(self):
        corr_scores = []
        for protein in tqdm(self.unique_proteins):
            query = get_entries_with_protein(protein, self.training_data)
            if len(query) == 1:
                corr_score = "NOT APPLICABLE"

            else:
                corr_score = get_corr_values(query)

            corr_scores.append(corr_score)

        return corr_scores

    def export_data(
            self,
            file_name: str,
            file_extension='csv',
            overwrite=False
    ) -> None:
        export_data_helper(self.unique_proteins_corr_data, file_name, file_extension, overwrite)


class InteractorJustifier:
    def __init__(self, training_data_path):
        self.training_data = pd.read_csv(training_data_path)
        self.unique_pairs = sorted(
            set(zip(self.training_data["UniProt_ID"], self.training_data["Mutation"]))
        )
        self.unique_pairs_corr_data = self.get_unique_pairs_correlation_data()

    def get_unique_pairs_correlation_data(self):
        corr_scores = self.get_corr_scores()
        unique_pairs_corr_data = pd.DataFrame(
            self.unique_pairs, columns=["PROTEIN", "MUTATION"]
        )
        unique_pairs_corr_data["PEARSON_CORR"] = corr_scores
        return unique_pairs_corr_data

    def get_corr_scores(self):
        corr_scores = []
        for pair in tqdm(self.unique_pairs):
            query = get_entries_with_pair(pair[0], pair[1], self.training_data)
            if is_same_interactor_pair(query):
                corr_score = "NOT APPLICABLE"
            else:
                corr_score = get_corr_values(query)

            corr_scores.append(corr_score)

        return corr_scores

    def export_data(
            self,
            file_name: str,
            file_extension='csv',
            overwrite=False
    ) -> None:

        export_data_helper(self.unique_pairs_corr_data, file_name, file_extension, overwrite)
