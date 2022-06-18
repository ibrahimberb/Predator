import pandas as pd
import numpy as np
from .displayers import display_labels

from ..mylogger import get_handler
import logging

from ..labels import ClassLabels

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.INFO)


def get_initial_columns_train(initial_columns_path):
    initial_features_train = pd.read_csv(initial_columns_path)["0"].to_list()
    triplet_columns = ["UniProt_ID", "Mutation", "Interactor_UniProt_ID"]
    # ['Mutation_Effect_Label'] + triplet columns + remaining feature list
    initial_columns_train_with_triplet = (
        [initial_features_train[0]] + triplet_columns + initial_features_train[1:]
    )

    return initial_columns_train_with_triplet


def get_initial_columns_tcga(initial_columns_path):
    initial_features_target = pd.read_csv(initial_columns_path)["0"].to_list()
    initial_features_target.remove("Mutation_Effect_Label")
    triplet_columns = ["UniProt_ID", "Mutation", "Interactor_UniProt_ID"]
    initial_columns_target = triplet_columns + initial_features_target

    return initial_columns_target


def reduce_columns_train(train_data, initial_columns_train_with_triplet):
    # Declare train data using selected column names
    train_data = train_data[initial_columns_train_with_triplet].copy()
    train_data.drop_duplicates(keep="first", inplace=True)
    log.debug(f"Size of dataframe: {train_data.shape}")
    log.debug("Dataframe head: {}".format(train_data.head()))
    assert train_data[train_data.duplicated()].empty

    return train_data


def reduce_columns_tcga(tcga_data, initial_columns_target):
    # Declare target_brca data using selected column names
    tcga_data = tcga_data[initial_columns_target].copy(deep=True)
    log.debug(f"Size of dataframe: {tcga_data.shape}")
    log.debug("Dataframe head: {}".format(tcga_data.head()))

    return tcga_data


def mutation_effect_label_binner(train_data, mutation_effect_label_column_name="Mutation_Effect_Label"):
    """
    Mutation Effect label binning is only applied to train_data.
    Apply Label binning.
        - Disruptive → ClassLabels.DISRUPTING
        - No effect + Increasing → ClassLabels.NONDISRUPTING
        - Decreasing → dropped
        - Causing → dropped
    """
    # Displaying possible label categories.
    if log.level == logging.DEBUG:
        display_labels(train_data)

    labels_to_bins = {
        "mutation disrupting(MI:0573)": ClassLabels.DISRUPTING,
        "mutation decreasing(MI:0119)": "IGNORED",
        "mutation disrupting strength(MI:1128)": ClassLabels.DISRUPTING,
        "mutation decreasing strength(MI:1133)": "IGNORED",
        "mutation with no effect(MI:2226)": ClassLabels.NONDISRUPTING,
        "disrupting": ClassLabels.DISRUPTING,
        "mutation increasing(MI:0382)": ClassLabels.NONDISRUPTING,
        "mutation increasing strength(MI:1132)": ClassLabels.NONDISRUPTING,
        "mutation decreasing rate(MI:1130)": "IGNORED",
        "mutation disrupting rate(MI:1129)": ClassLabels.DISRUPTING,
        "mutation causing(MI:2227)": "IGNORED",
        "mutation increasing rate(MI:1131)": ClassLabels.NONDISRUPTING,
    }

    replace_map = {mutation_effect_label_column_name: labels_to_bins}

    # Size of dataframe before binning.
    log.debug(f"Size of dataframe before binning: {train_data.shape}")

    # Modifications will be done on train_data_binned.
    train_data_binned = train_data.copy()

    # Replace the labels as described above.
    train_data_binned.replace(replace_map, inplace=True)

    # Drop the entries with "IGNORED": 'mutation cusing' in this case.
    train_data_binned = train_data_binned[
        train_data_binned[mutation_effect_label_column_name] != "IGNORED"
    ]

    # Reset index of the dataframe to avoid any possible errors
    train_data_binned.reset_index(drop=True, inplace=True)

    # Size of dataframe after binning.
    log.debug(f"Size of dataframe after binning: {train_data_binned.shape}")

    # First 5 rows of binned data.
    log.debug("Train Data Binned:")
    log.debug("Dataframe head: {}".format(train_data_binned.head()))

    # Confirming replacement of values are properly done. Mutation_Effect_Label only contains of 0 or 1.
    assert set(train_data_binned[mutation_effect_label_column_name].value_counts().index) == {0, 1}

    return train_data_binned


def type_coercion_data(data):
    # Some columns have been interpreted as object type, eventhough they are actually numeric.
    log.debug("{}".format(set(data.dtypes)))

    # These non-numeric interpereted columns will be coerced.  NaN  values will be converted to  0 .
    features = [
        column
        for column in data.columns
        if column not in ["UniProt_ID", "Mutation", "Interactor_UniProt_ID"]
    ]

    # Get column names where its type is *not* int or float, i.e. whose type is object.
    coerce_numeric_cols = set(
        [cname for cname in features if data[cname].dtype not in ["int64", "float64"]]
    )

    # Remove target variable from the list
    coerce_numeric_cols = coerce_numeric_cols - {
        "Mutation_Effect_Label",
        "UniProt_ID",
        "Mutation",
        "Interactor_UniProt_ID",
    }

    for cname in coerce_numeric_cols:
        data[cname] = pd.to_numeric(data[cname], errors="coerce")

    data = data.fillna(0)

    # Now all columns are interpreted as numeric type, except "UniProt_ID", "Mutation", "Interactor_UniProt_ID".
    log.debug("{}".format(set(data[features].dtypes)))
    assert set(data.dtypes[features].values) == {np.dtype("int64"), np.dtype("float64")}

    return data


def preprocess_train_data(train_data, initial_columns_path):
    initial_columns_train_with_triplet = get_initial_columns_train(initial_columns_path)
    train_data = reduce_columns_train(train_data, initial_columns_train_with_triplet)
    train_data_binned = mutation_effect_label_binner(train_data)
    train_data_processed = type_coercion_data(train_data_binned)
    return train_data_processed


def preprocess_tcga_data(tcga_data, initial_columns_path):
    initial_columns_tcga = get_initial_columns_tcga(initial_columns_path)
    tcga_data = reduce_columns_tcga(tcga_data, initial_columns_tcga)
    target_tcga_data = type_coercion_data(tcga_data)
    return target_tcga_data

