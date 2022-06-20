from datetime import datetime
import os.path as op
import numpy as np
from sklearn.base import BaseEstimator

from src.helpers.labels import ClassLabels


class BaselineProveanClassifier(BaseEstimator):
    """
    Prediction: "Deleterious" if the score is less than or equal to a predefined threshold.
    For our classification:
        if Provean <= thr:
            deleterious (disrupting)
        else:
            non deleterious (increasing or no effect)
    """
    PROVEAN_DEFAULT_THRESHOLD = -2.5
    FEATURE_NAME = "Provean_score"

    def __init__(self, feature_name=FEATURE_NAME, threshold=PROVEAN_DEFAULT_THRESHOLD):
        self.feature_name = feature_name
        self.threshold = threshold

    def fit(self, X, y=None):
        pass

    def predict(self, X):
        return np.array([self.apply_threshold(val) for val in X[self.feature_name]])

    def apply_threshold(self, val):
        if val <= self.threshold:
            return ClassLabels.DISRUPTING  # disrupting

        else:
            return ClassLabels.NONDISRUPTING  # increasing + no effect


def save_prediction_data(benchmark_dir, prediction_file_name, prediction_data):
    file_date = datetime.now().strftime("%Y-%m-%d")
    prediction_file_name = "{}_{}.csv".format(prediction_file_name, file_date)
    prediction_data.to_csv(op.join(benchmark_dir, prediction_file_name), index=False)
    print("Prediction data `{}`is exported.".format(op.join(benchmark_dir, prediction_file_name)))
