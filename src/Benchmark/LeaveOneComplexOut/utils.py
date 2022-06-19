import pandas as pd
from pandas import DataFrame
import numpy as np
import os.path as op
import shap
from sklearn.ensemble import RandomForestClassifier

from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import RepeatedStratifiedKFold

from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
)


from typing import List

PARAM_GRID_LEVEL_0 = {
    'max_depth': [2, 5, 10],
    'n_estimators': [10, 25, 50, 75, 100, 200, 400],
    'min_samples_split': [2, 5],
    'max_features': ['auto', 'sqrt', None],
    'class_weight': ['balanced', None]
}


class Hyperparameters:
    param_grid = PARAM_GRID_LEVEL_0
    n_iter = 10
    n_repeats_cv = 3
    n_jobs = -2
    verbose = 0


def get_scoring_metrics(y_true, y_pred):
    metrics_data = pd.DataFrame(
        [
            accuracy_score(y_true, y_pred),
            balanced_accuracy_score(y_true, y_pred),
            f1_score(y_true, y_pred),
            matthews_corrcoef(y_true, y_pred),
            precision_score(y_true, y_pred),
            recall_score(y_true, y_pred),
        ],
        columns=["SCORE"],
        index=["ACCURACY", "BALANCED_ACCURACY", "F1", "MATTEWS_COR", "PRECISION", "RECALL"])

    return metrics_data


def get_predictions(X_train, y_train, X_test, tuning=False):

    if tuning:
        randomized_search = get_randomized_search()
        randomized_search.fit(X_train, y_train)
        clf_tuned = randomized_search.best_estimator_
        predictions = clf_tuned.predict(X_test)

        return list(predictions)

    else:
        clf = get_default_classifier()
        clf.fit(X_train, y_train)
        predictions = clf.predict(X_test)

        return list(predictions)


def get_predictions_prob(X_train, y_train, X_test, tuning=False):

    if tuning:
        randomized_search = get_randomized_search()
        randomized_search.fit(X_train, y_train)
        clf_tuned = randomized_search.best_estimator_
        predictions_prob = clf_tuned.predict_proba(X_test)

        return list(predictions_prob)

    else:
        clf = get_default_classifier()
        clf.fit(X_train, y_train)
        predictions_prob = clf.predict_proba(X_test)

        return list(predictions_prob)


def check_train_protein_not_in_test_data(train_data, test_data):
    train_proteins = sorted(train_data["UniProt_ID"].unique())
    [test_protein] = test_data["UniProt_ID"].unique()
    assert test_protein not in train_proteins


def get_train_and_test_data(
        dataframes: List[DataFrame], index: int
):
    """
    Train data: all dataframes except the dataframe at given index.
    Test data: the dataframe at the given index.
    """
    train_dataframes = [dataframes[num] for num in range(len(dataframes)) if num != index]
    train_data = pd.concat(train_dataframes, ignore_index=True)
    test_data = dataframes[index]

    # Ensure test protein is not in train data.
    check_train_protein_not_in_test_data(train_data, test_data)

    return train_data, test_data


def get_default_classifier():
    return RandomForestClassifier(random_state=42, n_jobs=-1)


def read_features(file_path):
    with open(file_path) as fin:
        lines = fin.readlines()
        features = [line.strip() for line in lines]

    return features


def read_unique_proteins(file_path):
    with open(file_path) as fin:
        lines = fin.readlines()
        proteins = [line.strip() for line in lines]

    return proteins


def load_protein_dataframes(folder_path, order_list):
    protein_dataframes = []

    for item in order_list:
        hidden_entries_file_path = op.join(folder_path, f"hidden_entries_{item}.txt")
        data = pd.read_csv(hidden_entries_file_path, sep="\t")
        protein_dataframes.append(data)

    return protein_dataframes


def get_top_n_features(X_train, y_train, top_n):
    shap_values = get_shap_values(X_train, y_train)
    selected_features = get_selected_features(shap_values[1], X_train, top_n=top_n)
    print(f"{selected_features=}")

    return selected_features


def get_shap_values(X_train, y_train):
    """Fitted on train sets, not on all train data."""
    model = get_default_classifier()
    model.fit(X_train, y_train)
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_train, approximate=False, check_additivity=False)

    return shap_values


def get_selected_features(
        shap_values,
        features_data: DataFrame,
        top_n: int
) -> List[str]:
    column_list = features_data.columns
    feature_ratio = (np.abs(shap_values).sum(0) / np.abs(shap_values).sum()) * 100
    column_list = column_list[np.argsort(feature_ratio)[::-1]]
    column_list = column_list[:top_n]

    return list(column_list)


def get_randomized_search() -> RandomizedSearchCV:

    clf = get_default_classifier()

    randomized_search = RandomizedSearchCV(
        clf,
        Hyperparameters.param_grid,
        n_iter=Hyperparameters.n_iter,
        cv=RepeatedStratifiedKFold(
            n_splits=10, n_repeats=Hyperparameters.n_repeats_cv
        ),
        scoring='balanced_accuracy',
        return_train_score=True,
        n_jobs=Hyperparameters.n_jobs,
        verbose=Hyperparameters.verbose
    )

    return randomized_search
