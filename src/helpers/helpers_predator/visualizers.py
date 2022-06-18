from typing import List

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

from ..labels import ClassLabels

from ..mylogger import get_handler
import logging

handler_simple = get_handler('simple')
log_simple = logging.getLogger('Visualizers')
log_simple.handlers[:] = []
log_simple.addHandler(handler_simple)
log_simple.setLevel(logging.DEBUG)


def get_sampled_datasets_label_counts(sampled_train_data_list):
    sampled_train_data_to_label_counts = {}
    for i, sampled_train_data in enumerate(sampled_train_data_list):
        label_counts = sampled_train_data["Mutation_Effect_Label"].value_counts()
        sampled_train_data_to_label_counts["SAMPLED_TRAIN_DATA_" + str(i + 1)] = [
            label_counts.loc[0],
            label_counts.loc[1],
        ]
    sampled_datasets_label_counts = pd.DataFrame(sampled_train_data_to_label_counts).T
    sampled_datasets_label_counts.rename(
        {
            ClassLabels.DISRUPTING : 'Disrupting', 
            ClassLabels.NONDISRUPTING: 'NoEffect+Increasing'
        },
        axis="columns",
        inplace=True
    )

    sampled_datasets_label_counts.index.name = "SAMPLED_TRAIN_DATA"
    sampled_datasets_label_counts.reset_index(inplace=True)
    return sampled_datasets_label_counts


def visualize_sampled_train_datasets_label_counts(sampled_train_data_list, kind):

    sampled_datasets_label_counts = get_sampled_datasets_label_counts(
        sampled_train_data_list
    )

    experiment_statistics_data_melted = pd.melt(
        sampled_datasets_label_counts,
        id_vars=["SAMPLED_TRAIN_DATA"],
        value_vars=["NoEffect+Increasing", "Disrupting"],
        var_name="MUTATION_EFFECT",
        value_name="LABEL_COUNT",
    )

    medians = sorted(list(experiment_statistics_data_melted.groupby(['MUTATION_EFFECT'])['LABEL_COUNT'].median().values))
    medians = [int(median) for median in medians]
    print("MEDIANS:", medians)
    vertical_offset = experiment_statistics_data_melted['LABEL_COUNT'].median() * 0.05  # offset from median for display

    if kind in ["strip", "box"]:

        plt.figure(figsize=(3, 4))

        if kind == "strip":
            plot = sns.stripplot(
                x="MUTATION_EFFECT",
                y="LABEL_COUNT",
                data=experiment_statistics_data_melted,
                palette="ch:s=-.2,r=.6",
                jitter=True,
            )

        elif kind == "box":
            plot = sns.boxplot(
                x="MUTATION_EFFECT",
                y="LABEL_COUNT",
                data=experiment_statistics_data_melted,
                palette="ch:s=-.2,r=.6",
            )

        else:
            raise

        xticks = plot.get_xticks()
        plot.text(xticks[0], medians[xticks[0]] + vertical_offset, medians[xticks[0]],
                  horizontalalignment='left', size='small', color='k')
        plot.text(xticks[1], medians[xticks[1]] + vertical_offset, medians[xticks[1]],
                  horizontalalignment='right', size='small', color='k')

    elif kind == "bar":
        sampled_datasets_label_counts.plot(
            figsize=(25, 4), kind="bar", color=["#E3D9C1", "#27213F"],
            title='Label Counts per Experiment', xlabel='Experiment', ylabel='Counts',
            rot=0
        )

    else:
        log_simple.error(f"Parameter `kind` must be either `strip` or `bar`, not `{kind}`")


def visualize_accuracy_metrics(
    acc_scores: List[float],
    balan_acc_scores: List[float],
    kind='strip'
):
    sns.set_theme(style="ticks", palette="pastel", font_scale=1.15)
    data = pd.DataFrame({'ACCURACY': acc_scores,
                         'BALANCED_ACCURACY': balan_acc_scores})
    data_melted = pd.melt(data, var_name='METRIC', value_name='SCORE')
    # plt.figure(figsize=(4, 5))
    sns.catplot(x='METRIC', y='SCORE', data=data_melted,
                palette=sns.color_palette(['#265191', '#9F2945']),
                kind=kind)
    plt.show()


def visualize_distribution_top_n_features(shap_feature_selector, top_n):
    # plt.figure(figsize=(12, 12))  # poster purpose
    # sns.set(style='white', font_scale=2.5)  # poster purpose
    feature_to_counts = {}
    for feature in shap_feature_selector.n_features_to_aggregated_features[top_n]:
        count = (shap_feature_selector.aggregated_feature_selector
                 .n_features_to_selected_features_occurrences_counts[top_n])[feature]
        feature_to_counts[feature] = count

    counts_data = pd.DataFrame(feature_to_counts, index=['counts'])
    counts_data = pd.melt(counts_data, var_name='FEATURES', value_name='COUNTS')
    sns.barplot(x='FEATURES', y='COUNTS', color="#4C4C4C", data=counts_data)
    plt.title(f'Distribution of top-{top_n} features', fontweight="bold")
    plt.xlabel(None)
    plt.xticks(ha='right', rotation=45)

    sns.despine(right=True)
    # plt.tight_layout() # poster purpose
    # plt.savefig('foo.png')  # poster purpose
    plt.show()
