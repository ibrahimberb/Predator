from sklearn import metrics
import pandas as pd
import numpy as np


def get_fpr_tpr_ths(y_param, scores_param):
    """
    TODO: add positive_label parameter instead of hard coded value.
    Returns fpr, tpr, thresholds.
    Positive label is +.
    """

    y = np.array(y_param)
    scores = np.array(scores_param)
    fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label='+')

    return fpr, tpr, thresholds


def calculate_personalized_auc(data: pd.DataFrame):
    """
    Given personalized counts dataframe, calculate area under the curve (AUC) value
    of each columns with respect to target column. # TODO what is target column?? is having parameter really needed?
                                                   # it is CGC_STATUS, but having parameter may not be needed for now.
    Parameters
    ----------
    # TODO: add parameter doc string.


    Returns
    -------
    # TODO: ...
    auc_baseline, auc_ours, auc_ours_norm

    """

    # TODO: may be important: A temporary solution.
    #   Proteins were converted to gene names.
    #   however, not all of them were found correctly (about 15 of them), e.g. P123 → NA
    #   here, i am dealing with GENEs directly (without proteins), so is filled like this:
    #                    GENE  BASELINE  OUR_METHOD  OUR_METHOD_NORMALIZED CGC_STATUS
    #            15     GDF9         2           2                    1.0          -
    #            16    HSPA9         1           1                    1.0          -
    #            17  CSNK1A1         2           2                    1.0          -
    #   here ->  18       NA         0           0                    NaN          -
    #            19    LATS1         1           1                    1.0          +
    #   For now, I am just skipping those...
    data['OUR_METHOD_NORMALIZED'] = data['OUR_METHOD_NORMALIZED'].fillna(0)

    baseline_vs_cgc_status_data = data[["BASELINE", "CGC_STATUS"]].copy()
    ours_vs_cgc_status_data = data[["OUR_METHOD", "CGC_STATUS"]].copy()
    ours_normalized_vs_cgc_status_data = data[["OUR_METHOD_NORMALIZED", "CGC_STATUS"]].copy()

    fpr_baseline, tpr_baseline, _ = get_fpr_tpr_ths(np.array(baseline_vs_cgc_status_data["CGC_STATUS"]),
                                                    np.array(baseline_vs_cgc_status_data["BASELINE"]))

    fpr_ours, tpr_ours, _ = get_fpr_tpr_ths(np.array(ours_vs_cgc_status_data["CGC_STATUS"]),
                                            np.array(ours_vs_cgc_status_data["OUR_METHOD"]))

    fpr_ours_norm, tpr_ours_norm, _ = get_fpr_tpr_ths(np.array(ours_normalized_vs_cgc_status_data["CGC_STATUS"]),
                                                      np.array(ours_normalized_vs_cgc_status_data["OUR_METHOD_NORMALIZED"]))

    auc_baseline = metrics.auc(fpr_baseline, tpr_baseline)
    auc_ours = metrics.auc(fpr_ours, tpr_ours)
    auc_ours_norm = metrics.auc(fpr_ours_norm, tpr_ours_norm)

    return round(auc_baseline, 3), round(auc_ours, 3), round(auc_ours_norm, 3)
