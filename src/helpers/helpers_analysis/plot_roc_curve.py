from dataclasses import dataclass
from typing import Tuple

import numpy as np
from pandas import DataFrame
import pandas as pd

from sklearn import metrics
import matplotlib.pyplot as plt
import seaborn as sns

from .get_fpr_tpr_ths import get_fpr_tpr_ths

from ..mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)




class RocValues:
    def __init__(
            self,
            state_variables: np.ndarray,
            test_variables: np.ndarray
    ):
        """
        :param state_variables: E.g. GCG or CancerMine represented by '+'
        :param test_variables: E.g. "baseline" or  "our_method" column. values
        """
        fpr, tpr, ths = get_fpr_tpr_ths(
            np.array(state_variables), np.array(test_variables)
        )
        roc_auc = metrics.auc(fpr, tpr)
        self.fpr = fpr
        self.tpr = tpr
        self.ths = ths
        self.roc_auc = roc_auc


def plot_roc_curve_analysis(
    reference_data_name: str, cohort_specific,
    baseline_roc_values: RocValues,
    our_method_roc_values: RocValues,
    elaspic_cov_roc_values: RocValues
):

    plt.figure(figsize=(8, 6.5), dpi=600)
    sns.set_theme(style="white", palette="Set2", font_scale=2.00)  # Set2, Accent
    plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')

    if cohort_specific is None:
        reference_data_name = f"{reference_data_name.upper()}"
    else:
        reference_data_name = f"{reference_data_name.upper()}\ {cohort_specific.upper()}"

    # plt.title(f'Receiver Operating Characteristic (ROC)\n${reference_data_name}\ STATUS$ vs Various Columns')
    # TODO: line width: make it thicker
    plt.title(f'Receiver Operating Characteristic (ROC)\n${reference_data_name}$')  # STATUS
    plt.plot(baseline_roc_values.fpr, baseline_roc_values.tpr,
             label='baseline (%0.3f)' % baseline_roc_values.roc_auc)  # color="#d62728",
    plt.plot(our_method_roc_values.fpr, our_method_roc_values.tpr,
             label='our_method (%0.3f)' % our_method_roc_values.roc_auc)  # color="#bcbd22",
    plt.plot(elaspic_cov_roc_values.fpr, elaspic_cov_roc_values.tpr,
             label='elaspic_cov (%0.3f)' % elaspic_cov_roc_values.roc_auc)  # color="#17becf",

    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate\n(Sensitivity)')
    plt.xlabel('False Positive Rate\n(1-Specificity)')
    plt.show()

    log.debug("AUC BASELINE: {:.3f}".format(baseline_roc_values.roc_auc))
    log.debug("AUC OURS: {:.3f}".format(our_method_roc_values.roc_auc))
    log.debug("AUC ELASPIC COVERAGE: {:.3f}".format(elaspic_cov_roc_values.roc_auc))


def roc_curve_analysis(
        reference_data_name,
        preliminary_data,
        ref_gene_column,
        cohort_specific
):
    baseline_counts_vs_cgc_status_data = preliminary_data[["BASELINE", ref_gene_column]].copy()
    our_method_counts_vs_cgc_status_data = preliminary_data[["OUR_METHOD", ref_gene_column]].copy()
    elaspic_cov_vs_cgc_status_data = preliminary_data[["ELASPIC_COVERAGE", ref_gene_column]].copy()

    baseline_roc_values = RocValues(
        np.array(baseline_counts_vs_cgc_status_data[ref_gene_column]),
        np.array(baseline_counts_vs_cgc_status_data["BASELINE"])
    )

    our_method_roc_values = RocValues(
        np.array(our_method_counts_vs_cgc_status_data[ref_gene_column]),
        np.array(our_method_counts_vs_cgc_status_data["OUR_METHOD"])
    )

    elaspic_cov_roc_values = RocValues(
        np.array(elaspic_cov_vs_cgc_status_data[ref_gene_column]),
        np.array(elaspic_cov_vs_cgc_status_data["ELASPIC_COVERAGE"])
    )

    plot_roc_curve_analysis(
        reference_data_name=reference_data_name,
        cohort_specific=cohort_specific,
        baseline_roc_values=baseline_roc_values,
        our_method_roc_values=our_method_roc_values,
        elaspic_cov_roc_values=elaspic_cov_roc_values
    )

    auc_scores = {
        "BASELINE": "%0.3f" % baseline_roc_values.roc_auc,
        "OURS": "%0.3f" % our_method_roc_values.roc_auc,
        "ELASPIC_COV": "%0.3f" % elaspic_cov_roc_values.roc_auc
    }

    return auc_scores

    # return baseline_roc_values, our_method_roc_values, plot_roc_curve_analysis # ?????

    # # displayers
    # display(pd.DataFrame({
    #     "roc_auc_baseline": [roc_auc_baseline],
    #     "roc_auc_our_method": [roc_auc_our_method],
    #     "roc_auc_elaspic_cov": [roc_auc_elaspic_cov],
    # }, index=['AUC']).T)
    #
    # display(pd.DataFrame({
    #     "roc_auc_baseline_brca": [roc_auc_baseline_tcga],
    #     "roc_auc_our_method_brca": [roc_auc_our_method_tcga],
    #     "roc_auc_elaspic_cov_brca": [roc_auc_elaspic_cov_tcga],
    # }, index=['AUC']).T)
    #
    # display(pd.DataFrame({
    #     "roc_auc_baseline_bnz": [roc_auc_baseline_bnz],
    #     "roc_auc_our_method_bnz": [roc_auc_our_method_bnz],
    #     "roc_auc_elaspic_cov_bnz": [roc_auc_elaspic_cov_bnz],
    # }, index=['AUC']).T)
    #
    # display(pd.DataFrame({
    #     "roc_auc_baseline_bnz_brca": [roc_auc_baseline_bnz_brca],
    #     "roc_auc_our_method_bnz_brca": [roc_auc_our_method_bnz_brca],
    #     "roc_auc_elaspic_cov_bnz_brca": [roc_auc_elaspic_cov_bnz_brca],
    # }, index=['AUC']).T)