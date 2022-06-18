from sklearn import metrics
import numpy as np


def get_fpr_tpr_ths(y_param, scores_param):
    """
    Returns fpr, tpr, thresholds.
    Positive label is +.
    """

    y = np.array(y_param)
    scores = np.array(scores_param)
    fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label='+')

    return fpr, tpr, thresholds
