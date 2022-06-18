from dataclasses import dataclass
from typing import List, Dict

from .data_materials import DataMaterials
from .machine_learning_utils import get_default_classifier

from tqdm.notebook import tqdm
from collections import defaultdict

from pandas import DataFrame

import numpy as np
import pandas as pd

import shap

from ..mylogger import get_handler
import logging

handler = get_handler()
handler_simple = get_handler('simple')

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)

log_simple = logging.getLogger('Feature_selection')
log_simple.handlers[:] = []
log_simple.addHandler(handler_simple)
log_simple.setLevel(logging.DEBUG)


@dataclass
class FeatureAggregationMethods:
    OCCURRENCE = 'occurrence'
    CONCORDANCE = 'concordance'


class ShapFeatureSelector:
    def __init__(
            self,
            n_experiment,
            shap_top_ns: List[int]
    ):
        log.debug("Initializing ShapFeatureSelector ..")
        self.shap_top_ns = shap_top_ns
        self.n_experiment = n_experiment
        shap.initjs()

        self.shap_values_train_list = None

        self.aggregated_feature_selector = None

        self.n_features_to_selected_features_list = defaultdict(list)  # Dict[int, List[List[str]]]
        self.n_features_to_aggregated_features = None

    def load_shap_values(self, data_materials: DataMaterials):
        """Fitted on train sets, not on all train data."""
        log.debug("Loading ShapFeatureSelector ..")
        shap_values_train_list = []
        for i in tqdm(range(self.n_experiment)):
            model = get_default_classifier()
            model.fit(data_materials["Xs_train"][i], data_materials["ys_train"][i])
            explainer = shap.TreeExplainer(model)
            shap_values_train_list.append(explainer.shap_values(data_materials["Xs_train"][i],
                                                                approximate=False, check_additivity=False))

        self.shap_values_train_list = shap_values_train_list

    @staticmethod
    def get_selected_features_single_data(shap_values, features_data: DataFrame, top_n) -> List[str]:
        column_list = features_data.columns
        feature_ratio = (np.abs(shap_values).sum(0) / np.abs(shap_values).sum()) * 100
        column_list = column_list[np.argsort(feature_ratio)[::-1]]
        column_list = column_list[:top_n]

        return list(column_list)

    def get_selected_features(self, data_materials: DataMaterials):
        log_simple.debug(" === SELECTED FEATURES === ")
        for top_n in self.shap_top_ns:
            log_simple.debug(f" --- SHAP TOP {top_n} ---")
            for exp in range(self.n_experiment):
                shap_values = self.shap_values_train_list[exp][1]
                features_data = data_materials["Xs_train"][exp]
                selected_features = self.get_selected_features_single_data(shap_values, features_data, top_n=top_n)
                log_simple.debug(f"Experiment {exp + 1}")
                print(f"{selected_features}\n")
                self.n_features_to_selected_features_list[top_n].append(selected_features)

    def shap_aggregate_selected_features(self, method):
        self.aggregated_feature_selector = AggregatedFeatureSelector(self.n_features_to_selected_features_list, method)
        self.aggregated_feature_selector.aggregate()
        self.n_features_to_aggregated_features = self.aggregated_feature_selector.n_features_to_aggregated_features

    def display_rankings(self, top_n, extract):
        log.debug("Displaying rankings ..")
        data = pd.DataFrame(self.n_features_to_selected_features_list[top_n])
        filename = f'rankings_{top_n}.csv'
        if extract:
            log.debug("Extracting rankings ..")
            data.to_csv(filename, index=False)
            log.info(f"Rankings are extracted to file {filename}.")

        return data.head()

    def select_features(self, data_materials: DataMaterials):
        self.load_shap_values(data_materials)


class AggregatedFeatureSelector:
    def __init__(self, n_features_to_selected_features_list: Dict[int, List[List[str]]], aggregation_method=None):
        self.aggregation_method = aggregation_method
        self.n_features_to_selected_features_list = n_features_to_selected_features_list
        self.n_features_to_aggregated_features = None
        self.n_features_to_selected_features_occurrences_counts = defaultdict(dict)

    def aggregate(self):
        if self.aggregation_method == FeatureAggregationMethods.OCCURRENCE:
            self.n_features_to_aggregated_features = self.aggregate_occurrence()

        else:
            log.critical("Aggregation method not recognized.")
            raise ValueError(f"Invalid Aggregation method: {self.aggregation_method}")

    def aggregate_occurrence(self) -> Dict[int, list]:
        top_n_to_frequently_occurred_features = {}
        for top_n, features_list in self.n_features_to_selected_features_list.items():
            selected_features_to_occurrences = defaultdict(int)
            for features in features_list:
                for feature in features:
                    selected_features_to_occurrences[feature] += 1

            # Sorting the counts dictionary by value in descending order.
            frequently_occurred_features = sorted(selected_features_to_occurrences,
                                                  key=selected_features_to_occurrences.get, reverse=True)[:top_n]
            top_n_to_frequently_occurred_features[top_n] = frequently_occurred_features
            self.n_features_to_selected_features_occurrences_counts[top_n] = selected_features_to_occurrences

        return top_n_to_frequently_occurred_features
