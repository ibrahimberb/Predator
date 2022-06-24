from typing import List

from pandas import DataFrame
import pandas as pd
import numpy as np

from ..mylogger import get_handler
import logging

from tqdm.notebook import tqdm

from ..labels import ClassLabels

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


def get_predictive_columns_removed_data(data: DataFrame) -> DataFrame:
    """
    Remove predictive columns (i.e. feature columns) from given dataset, leaving them
    with predictions and triplets.

    Parameters
    ----------
        data : <DataFrame>
            The input dataframe whose predictive columns to be removed.

    Returns
    -------
        features_removed_data : <DataFrame>
            The dataframe containing prediction values along with triplet information.
            I.e. ["Prediction", "UniProt_ID", "Mutation", "Interactor_UniProt_ID"]
    """

    features_removed_data = data[["Prediction", "UniProt_ID", "Mutation", "Interactor_UniProt_ID"]].copy()
    return features_removed_data


def max_votes(x):
    """
    Return the maximum occurrence of predicted class.

    Notes
    -----
        If number of class ClassLabels.DISRUPTING prediction is equal to number of class ClassLabels.NONDISRUPTING predictions, 
        NO_VOTE will be returned.
        E.g.
            Num_preds_ClassLabels.DISRUPTING = 25,
            Num_preds_ClassLabels.NONDISRUPTING = 25,
            Num_preds_NO_VOTE = 0,
            returned vote : "NO_VOTE".

    """
    if x[f'Num_preds_{ClassLabels.DISRUPTING}'] > x[f'Num_preds_{ClassLabels.NONDISRUPTING}'] and x[f'Num_preds_{ClassLabels.DISRUPTING}'] > x['Num_preds_NO_VOTE']:
        return ClassLabels.DISRUPTING
    elif x[f'Num_preds_{ClassLabels.NONDISRUPTING}'] > x[f'Num_preds_{ClassLabels.DISRUPTING}'] and x[f'Num_preds_{ClassLabels.NONDISRUPTING}'] > x['Num_preds_NO_VOTE']:
        return ClassLabels.NONDISRUPTING
    else:
        return 'NO_VOTE'


def convert_primary_isomer(column_name: str, data: DataFrame) -> DataFrame:
    """
    Converts proteins into primary form representation (dash-free from) in given column name of the given dataframe.
    E.g.
        P16473-2 → P16473

    Parameters
    ----------
        column_name : <string>
            Name of the column where protein is stored.

        data : <DataFrame>
            The dataframe whose proteins will be processed in `column_name` column.

    Returns
    -------
        data : <DataFrame>
            Processed version of input dataframe.

    """

    # Protein names will be converted dashed-free version, if they contain.
    data[column_name] = data[column_name].apply(lambda x: x.split('-')[0])

    return data


def get_prediction_entry(protein, mutation, interactor, data):
    predictions_array = data[
        (data['UniProt_ID'] == protein) &
        (data['Mutation'] == mutation) &
        (data['Interactor_UniProt_ID'] == interactor)
        ]['Prediction'].values

    # The entry such that it is predicted both `1` and `0` for that triplet, hence dropped.
    if predictions_array.size == 0:
        return 'NO_VOTE'

    # Prediction array contains one element, and return that element.
    elif len(predictions_array) == 1:
        [prediction] = predictions_array  # extracting single value from length-1 list.
        return prediction

    else:
        raise ValueError('There should be one entry, thus one prediction value, '
                         'but contains {} elements.'.format(len(predictions_array)))


def get_prediction_entry_prob(row, data):
    """
    Returns the probability of being class X. # TODO: Which class?
    :param row:
    :param data:
    :return:
    """
    prediction_probs_array = data[
        (data['UniProt_ID'] == row['UniProt_ID']) &
        (data['Mutation'] == row['Mutation']) &
        (data['Interactor_UniProt_ID'] == row['Interactor_UniProt_ID'])
    ]['Prediction'].values

    # log.info(f"prediction_probs_array: \n{prediction_probs_array}")

    # The entry such that it is predicted both higher and lower than 0.50 for that triplet, hence dropped.
    if prediction_probs_array.size == 0:
        return 'NO_VOTE'

    # Prediction array contains one element, and return that element.
    elif len(prediction_probs_array) == 1:
        [prediction] = prediction_probs_array  # extracting probability of class X # TODO: Which class?
        return prediction

    # Prediction array contains more elements, but they should be on the same probability intervals
    # ie. all being lower than 0.50 or all being above 0.50.
    elif are_all_valid_predictions(prediction_probs_array):
        # log.warning(f"are_all_valid_predictions :\n{prediction_probs_array}")
        return np.mean(prediction_probs_array)

    else:
        log.critical(f"prediction_probs_array: {prediction_probs_array}")
        raise ValueError('** prob ** There should be one entry, thus one prediction value, '
                         'but contains {} elements.'.format(len(prediction_probs_array)))


def get_triplet_columns(data: DataFrame) -> DataFrame:
    """
    Remove all columns except (protein, mutation, interactor) triplets.

    Parameters
    ----------
        data : <DataFrame>
            The input dataframe whose non-triplet columns to be removed.

    Returns
    -------
        triplet_data : <DataFrame>
            The dataframe containing triplet information with following columns:
            ["UniProt_ID", "Mutation", "Interactor_UniProt_ID"]
    """

    triplet_data = data[["UniProt_ID", "Mutation", "Interactor_UniProt_ID"]].copy()
    return triplet_data


def remove_triplet_columns(data: DataFrame) -> DataFrame:
    """
    Remove all columns except (protein, mutation, interactor) triplets.

    Parameters
    ----------
        data : <DataFrame>
            The input dataframe whose non-triplet columns to be removed.

    Returns
    -------
        triplet_data : <DataFrame>
            The dataframe containing triplet information with following columns:
            ["UniProt_ID", "Mutation", "Interactor_UniProt_ID"]
    """

    triplet_removed_data = data.drop(["UniProt_ID", "Mutation", "Interactor_UniProt_ID"], axis='columns').copy()
    return triplet_removed_data


def drop_invalid_predicted_entries(data: DataFrame):
    """
    Prediction data contains entries which for the same (PROTEIN, MUTATION, INTERACTOR), the predicted
    class is both ClassLabels.DISRUPTING and ClassLabels.NONDISRUPTING. Find such instances, and drop them.

    Parameters
    ----------
        data : <DataFrame>
            The dataframe whose invalid predicted entries will be dropped.

    Returns
    -------
        data : <DataFrame>
            Processed version of input dataframe.

        removed_entries_data : <DataFrame>
            A dataframe which contains removed entires.
    """

    entries = []
    entries_ix = []

    # For each (PROTEIN, MUTATION, INTERACTOR), capture the predicted class numbers.
    # If they are not all the same, then that (PROTEIN, MUTATION, INTERACTOR) row will be dropped.
    for index, row in data.iterrows():
        # Predicted class number(s) for current (PROTEIN, MUTATION, INTERACTOR) triplet.
        # Ideally, should be all the same. If not, then it will contain two class names.
        search_data_predictions = data[(data["UniProt_ID"] == row["UniProt_ID"]) &
                                       (data["Mutation"] == row["Mutation"]) &
                                       (data["Interactor_UniProt_ID"] == row["Interactor_UniProt_ID"])][
            "Prediction"].unique()

        # If search_data_predictions contains class-0 and class-1 together, then it is an invalid predicted entries.
        if len(search_data_predictions) > 1:
            entries.append((row["Prediction"], row["UniProt_ID"], row["Mutation"], row["Interactor_UniProt_ID"]))
            entries_ix.append(index)

    data_dropped, removed_entries_data = get_valid_and_invalid_entries_datasets(data, entries, entries_ix)

    return data_dropped, removed_entries_data


def are_all_valid_predictions(prediction_probs):
    return all(prediction_probs > 0.50) or all(prediction_probs <= 0.50)


def drop_invalid_predicted_probs_entries(data: DataFrame):

    entries = []
    entries_ix = []

    # For each (PROTEIN, MUTATION, INTERACTOR), capture the predicted class probabilities.
    # If all probabilities for particular triplet are not on the same probability portion,
    # (being above 0.50 or below), then that (PROTEIN, MUTATION, INTERACTOR) row will be dropped.
    for index, row in data.iterrows():
        # Predicted class number(s) for current (PROTEIN, MUTATION, INTERACTOR) triplet.
        # Ideally, should be all the same. If not, then it will contain two class names.
        search_data_predictions_probs = data[
            (data["UniProt_ID"] == row["UniProt_ID"]) &
            (data["Mutation"] == row["Mutation"]) &
            (data["Interactor_UniProt_ID"] == row["Interactor_UniProt_ID"])
            ]["Prediction"]

        # If search_data_predictions_probs contains probabilities values of greater than 0.50 and lower than 50
        # at the same time, then it is an invalid predicted entries.
        if not are_all_valid_predictions(search_data_predictions_probs):
            entries.append((row["Prediction"], row["UniProt_ID"], row["Mutation"], row["Interactor_UniProt_ID"]))
            entries_ix.append(index)

    data_dropped, removed_entries_data = get_valid_and_invalid_entries_datasets(data, entries, entries_ix)

    return data_dropped, removed_entries_data


def get_valid_and_invalid_entries_datasets(data, entries, entries_ix):
    removed_entries_data = pd.DataFrame(
        entries, columns=["PREDICTION", "PROTEIN", "MUTATION", "INTERACTOR"]
    )

    log.debug('Removed entries first five rows (of {}): \n{}'.format(
        removed_entries_data.shape[0], removed_entries_data.head())
    )

    # Drop invalid predicted entries based on their index.
    data_dropped = data.drop(entries_ix, axis='index')

    # Reset index of the dataframe to avoid any possible errors.
    data_dropped.reset_index(drop=True, inplace=True)

    return data_dropped, removed_entries_data


def add_votes(
        data: DataFrame,
        final_prediction_datasets: List[DataFrame]
) -> DataFrame:
    """
    Add the votes from final prediction datasets to given data.
    For each entry in data to become ensambled prediction data, the corresponding predictions will be looked up
    and placed at that entry. If corresponding prediction does not exists (indicating that prediction removed from
    that prediction data because it was an invalid prediction), "NO_VOTE" label will be assigned.

    Data will have the following form:
    # todo, table missing 3 columns, Num_preds_0, Num_preds_1, Num_preds_NO_VOTE ???????????????
        +---------+----------+------------+--------------+--------------+-----+--------------+
        | Protein | Mutation | Interactor | Prediction 1 | Prediction 2 | ... | Prediction N |
        +---------+----------+------------+--------------+--------------+-----+--------------+
        |         |          |            | 1            | 0            | ... | 1            |
        +---------+----------+------------+--------------+--------------+-----+--------------+
        |         |          |            | 0            | 0            | ... | NO_VOTE      |
        +---------+----------+------------+--------------+--------------+-----+--------------+
        |         |          |            | 1            | 1            | ... | 1            |
        +---------+----------+------------+--------------+--------------+-----+--------------+

    Parameters
    ----------
        data : <DataFrame>
            The data to be ensambled prediction data, which final predictions will be added to.

        final_prediction_datasets : <List[DataFrame]>
            A list of datasets containing predictions.

    Returns
    -------
        data : <DataFrame>
            Corresponding predictions added version of input data.
    """

    final_votes = []
    for index, row in tqdm(data.iterrows(), total=len(data)):
        protein, mutation, interactor = row['UniProt_ID'], row['Mutation'], row['Interactor_UniProt_ID']

        votes = []
        for final_prediction_data in final_prediction_datasets:
            prediction = get_prediction_entry(protein, mutation, interactor, final_prediction_data)
            votes.append(prediction)

        final_votes.append((votes.count(0), votes.count(1), votes.count('NO_VOTE')))

    data[f'Num_preds_{ClassLabels.DISRUPTING}'] = [item for item, _, _ in final_votes]
    data[f'Num_preds_{ClassLabels.NONDISRUPTING}'] = [item for _, item, _ in final_votes]
    data['Num_preds_NO_VOTE'] = [item for _, _, item in final_votes]

    return data


def take_avg(x):
    """
    Returns the average of probabilities.
    If number of NO_VOTE are more than half, returns NO_VOTE.
    """
    n_no_votes = len([e for e in x if e == 'NO_VOTE'])
    if n_no_votes >= (len(x) / 2):
        return "NO_VOTE"
    else:
        return round(np.mean([e for e in x if e != 'NO_VOTE']), 5)


def take_median(x):
    """
    Returns the median of probabilities.
    If number of NO_VOTE are more than half, returns NO_VOTE.
    """
    n_no_votes = len([e for e in x if e == 'NO_VOTE'])
    if n_no_votes >= (len(x) / 2):
        return "NO_VOTE"
    else:
        return round(np.median([e for e in x if e != 'NO_VOTE']), 5)


def add_voted_probs(
        data: DataFrame,
        final_prediction_datasets: List[DataFrame]
) -> DataFrame:
    """
    Add the vote probabilities from final prediction datasets to given data.
    For each entry in data to become ensemble prediction data, the corresponding prediction probabilities will be
    looked up, averaged, and placed at that entry. If corresponding prediction probability does not exists
    (indicating that prediction removed from that prediction data because it was an invalid prediction),
    "NO_VOTE" label will be assigned.

    Data will have the following form:
        +---------+----------+------------+--------------+--------------+-----+--------------+
        | Protein | Mutation | Interactor | Prediction 1 | Prediction 2 | ... | Prediction N |
        +---------+----------+------------+--------------+--------------+-----+--------------+
        |         |          |            | 0.65         | 0.35         | ... | 0.90         |
        +---------+----------+------------+--------------+--------------+-----+--------------+
        |         |          |            | 0.20         | 0.40         | ... | NO_VOTE      |
        +---------+----------+------------+--------------+--------------+-----+--------------+
        |         |          |            | 0.45         | 0.80         | ... | 0.90         |
        +---------+----------+------------+--------------+--------------+-----+--------------+

    Notes
    -----
        In taking the average, we exclude NO_VOTE. For example, the calculation of a sample entry
        will be as follows when number of experiment is 5:
            Pred_1      Pred_2      Pred_3      Pred_4      Pred_5
              0.20        0.40     NO_VOTE        0.30        0.80

            Avg = (0.20 + 0.40 + 0.30 + 0.80) / 4
                = 0.425

    Parameters
    ----------
        data : <DataFrame>
            The data to be ensemble prediction data, which final prediction probabilities will be added to.

        final_prediction_datasets : <List[DataFrame]>
            A list of datasets containing prediction probabilities.

    Returns
    -------
        data : <DataFrame>
            Corresponding predictions probabilities added version of input data.
    """

    for i, final_prediction_data in tqdm(
            enumerate(final_prediction_datasets), total=len(final_prediction_datasets)
    ):
        data[f'Trial {i}'] = data.apply(get_prediction_entry_prob, axis=1, args=(final_prediction_data,))
        # x = my_series.apply(my_function, args=(arg1,))

    return data
