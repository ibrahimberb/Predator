from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import cross_val_score

from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import make_scorer

import pandas as pd
import numpy as np
from tqdm.notebook import tqdm

from .machine_learning_utils import get_default_classifier, evaluate_valid

import seaborn as sns
from matplotlib.ticker import MultipleLocator
from matplotlib import pyplot as plt

from IPython.display import display

from ..mylogger import get_handler
import logging

handler = get_handler()
handler_simple = get_handler("simple")

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)

log_simple = logging.getLogger(__name__ + '1')
log_simple.handlers[:] = []
log_simple.addHandler(handler_simple)
log_simple.setLevel(logging.DEBUG)

mcc = make_scorer(matthews_corrcoef)

EVAL_METRICS_DICT = {
    "F1": "f1",
    "BALANCED_ACCURACY": "balanced_accuracy",
    "ACCURACY": "accuracy",
    "F1_MACRO": "f1_macro",
    "F1_MICRO": "f1_micro",
    "PRECISION": "precision",
    "RECALL": "recall",
    "ROC_AUC": "roc_auc",
    "PRECISION_MACRO": "precision_macro",
    "PRECISION_MICRO": "precision_micro",
    # Custom metrics
    "MCC": mcc,
}

# An option is to use `cross_validate` function. It takes multiple scoring metrics.
# https://stackoverflow.com/questions/35876508/evaluate-multiple-scores-on-sklearn-cross-val-score/43978914


class EvaluationMetrics:
    def __init__(self, n_experiment, data_materials, shap_feature_selector):
        log.info("Initializing EvaluationMetrics..")
        self.n_experiment = n_experiment
        self.shap_feature_selector = shap_feature_selector
        self.data_materials = data_materials
        self.n_repeats = None
        self.Xs_benchmark_feature_names_to_dataframes_list = None
        self.benchmark_columns = None
        self.scoring_metrics_list = None
        self.scoring_metrics_data = None
        self.scoring_metrics_data_melted = None

        self.initialize_benchmark_dataframes()

    def initialize_benchmark_dataframes(self):
        # We benchmark dataframes with train sets.
        log.info("Initialize_benchmark_dataframes ..")
        self.Xs_benchmark_feature_names_to_dataframes_list = []
        self.benchmark_columns = ['Provean', r'$\Delta\Delta$G', 'All Columns']
        for n in self.shap_feature_selector.shap_top_ns:
            self.benchmark_columns.append(f"SHAP Columns ({n})")

        for exp in range(self.n_experiment):
            Xs_benchmark_feature_names_to_dataframes = {
                f"X_train_exp_{exp}_provean": self.data_materials["Xs_train_provean"][exp],
                f"X_train_exp_{exp}_ddG": self.data_materials["Xs_train_ddG"][exp],
                f"X_train_exp_{exp}": self.data_materials["Xs_train"][exp]
            }
            for n in self.shap_feature_selector.shap_top_ns:
                Xs_benchmark_feature_names_to_dataframes.update({
                    f"X_train_exp_{exp}_shap_{n}": self.data_materials[f"Xs_train_shap_{n}"][exp]
                })

            self.Xs_benchmark_feature_names_to_dataframes_list.append(Xs_benchmark_feature_names_to_dataframes)

    def run_eval_metrics(self, n_repeats, n_jobs, verbose=False):
        log.info("Running evaluation metrics ..")
        self.n_repeats = n_repeats
        scoring_metrics_list = []
        for exp in tqdm(range(self.n_experiment)):
            scoring_metrics = evaluate_metrics(
                self.Xs_benchmark_feature_names_to_dataframes_list[exp],
                self.data_materials["ys_train"][exp],
                n_repeats=n_repeats,
                n_jobs=n_jobs,
                verbose=verbose
            )
            scoring_metrics_list.append(scoring_metrics)
        self.scoring_metrics_list = scoring_metrics_list
        self.prepare_scoring_metrics_data()

    def prepare_scoring_metrics_data_melted(self):
        # Melting the dataframe.
        scoring_metrics_entries = []
        for i in range(self.n_experiment):
            scoring = self.scoring_metrics_list[i].copy()
            scoring.index.name = 'X_NAME'
            scoring.insert(0, 'FEATURES', self.benchmark_columns)
            scoring.insert(0, 'EXPERIMENT_NO', i)
            scoring.reset_index(inplace=True)
            scoring = scoring.melt(id_vars=["X_NAME", "FEATURES", "EXPERIMENT_NO"],
                                   var_name="METRIC",
                                   value_name="SCORE")
            scoring_metrics_entries.append(scoring)

        scoring_metrics_data_melted = pd.concat(scoring_metrics_entries, axis='rows')
        print(scoring_metrics_data_melted.shape)
        scoring_metrics_data_melted.head()
        self.scoring_metrics_data_melted = scoring_metrics_data_melted

    def prepare_scoring_metrics_data(self):
        self.prepare_scoring_metrics_data_melted()
        dataframes = []
        for column in self.benchmark_columns:
            dataframes.append(
                (self.scoring_metrics_data_melted[self.scoring_metrics_data_melted['FEATURES'] == column]
                 .groupby(by=['METRIC'])
                 .mean()
                 .drop(['EXPERIMENT_NO'], axis='columns')
                 .rename(columns={'SCORE': column}))
            )
        self.scoring_metrics_data = pd.concat(dataframes, axis='columns')

    def plot_performance_comparison_results(self, metrics=None):
        # sns.set_theme(style="ticks", palette="pastel", font_scale=1.5)  # TODO: POSTER,uncommendLATER
        sns.set_theme(style="ticks", palette="pastel", font_scale=1.65)
        # TODO: [later] plot size adjusting itself depending on input ↓
        if metrics is None:
            plt.figure(figsize=(30 + len(self.benchmark_columns), 7))
        else:
            plt.figure(figsize=(3 * len(metrics), 7))
# title_string_1 = fr"Performance\ Comparison\ of\ Selected\ Features\ vs.\ All\ Features"
        # title_string_2 = fr"CV = 10, CV\_repeat = {self.n_repeats}, Experiment\_repeat = {self.n_experiment}"
        title_string_1 = fr"Performance\ Comparison\ of\ Selected\ Features"
        title_string_2 = ""
        plt.title(f"$\mathbf{{{title_string_1}}}$ \n $\mathbf{{{title_string_2}}}$", fontsize=24, fontweight='bold')
        plt.ylabel('Metrics', fontsize=24, fontweight='bold')
        plt.xlabel('Scores', fontsize=24, fontweight='bold')
        plt.axhline(y=0.5, color='k', linestyle='--', alpha=0.8, lw=0.5)
        # noinspection SpellCheckingInspection

        if metrics is not None:
            scoring_metrics_data_melted_less_metrics = self.scoring_metrics_data_melted[
                self.scoring_metrics_data_melted['METRIC'].isin(metrics)
            ]
            scoring_metrics_data_plot = scoring_metrics_data_melted_less_metrics

        else:
            scoring_metrics_data_plot = self.scoring_metrics_data_melted

        # ax = sns.boxplot(x='METRIC', y='SCORE', hue='FEATURES', data=scoring_metrics_data_melted_less_metrics,
        #                  palette='Pastel1')  # bone, vlag, cividis, #03012d, light:#444452

        ax = sns.boxplot(x='METRIC', y='SCORE', hue='FEATURES', data=scoring_metrics_data_plot,
                         palette='Pastel1')  # bone, vlag, cividis, #03012d, light:#444452
        ax.xaxis.set_minor_locator(MultipleLocator(0.5))
        ax.xaxis.grid(True, which='minor', color='#ababab', lw=1)
        # legend = plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), edgecolor="black") ############
        # legend.get_frame().set_alpha(None)
        ###################
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), framealpha=0)
        # plt.tight_layout()  # poster purpose
        # plt.savefig('foo2.png')  # poster purpose
        plt.show()
        # sns.despine(offset=10, trim=True)


class EvaluationValid:
    def __init__(
            self,
            n_experiment,
            data_materials,
            default_models=None,
            tuned_models=None,
    ):
        self.n_experiment = n_experiment
        self.data_materials = data_materials
        self.default_models = default_models
        self.tuned_models = tuned_models
        self.qualified_models = None
        self.scores = {}
        self.comparison_data = None

    def get_train_valid_splits(self, exp, determined_features):
        if determined_features is None:
            X_train = self.data_materials["Xs_train"][exp]
            X_valid = self.data_materials["Xs_valid"][exp]

        else:
            X_train = self.data_materials["Xs_train"][exp][determined_features]
            X_valid = self.data_materials["Xs_valid"][exp][determined_features]

        y_train = self.data_materials["ys_train"][exp]
        y_valid = self.data_materials["ys_valid"][exp]

        return X_train, y_train, X_valid, y_valid

    def run_evaluation(
            self,
            classifiers,
            show_confusion_matrix=False,
            determined_features=None
    ):
        acc_scores = []
        balan_acc_scores = []
        for exp in tqdm(range(self.n_experiment)):
            print("-------- EXPERIMENT: {:>2} --------".format(exp + 1))
            X_train, y_train, X_valid, y_valid = self.get_train_valid_splits(exp, determined_features)
            log_simple.debug(f"X_train.shape={X_train.shape}, y_train.shape={y_train.shape}, "
                             f"X_valid.shape={X_valid.shape}, y_valid.shape={y_valid.shape}")
            log_simple.debug(f"Classifier: {classifiers[exp]}")
            acc_score, balan_acc_score = evaluate_valid(
                classifiers[exp], X_train, y_train, X_valid, y_valid, show_confusion_matrix
            )
            acc_scores.append(acc_score)
            balan_acc_scores.append(balan_acc_score)
            print("================================")

        acc_scores_avg = np.array(acc_scores).mean()
        balan_acc_scores_avg = np.array(balan_acc_scores).mean()
        scores_dict = {
            "acc_scores": acc_scores,
            "balan_acc_scores": balan_acc_scores,
            "acc_scores_avg": acc_scores_avg,
            "balan_acc_scores_avg": balan_acc_scores_avg
        }
        return scores_dict

    def evaluate(self, models_type, show_confusion_matrix, determined_features):
        log.debug("Training on train set and measuring performance by predicting on validation set.")
        if models_type == "default":
            log.debug("Evaluating with default models.")
            scores_dict = self.run_evaluation(
                self.default_models, show_confusion_matrix
            )
            self.scores["initial_scoring"] = scores_dict

        elif models_type == "feature_selected":
            # Feature selected models are essentially default models
            # but we run then with determined features.
            log.debug("Evaluating with default models using determined features.")
            log.debug(f"Determined features: \n{determined_features}")
            scores_dict = self.run_evaluation(
                self.default_models, show_confusion_matrix, determined_features
            )
            self.scores["feature_selected_scoring"] = scores_dict

        elif models_type == "tuned":
            log.debug("Evaluating with tuned models.")
            log.debug(f"Determined features: \n{determined_features}")
            scores_dict = self.run_evaluation(
                self.tuned_models, show_confusion_matrix, determined_features
            )
            self.scores["finalized_scoring"] = scores_dict

        else:
            raise ValueError("Invalid arg for `models_type`.")

    def compare_models(self, kind):
        """
        Comparison of following models:
            1. Default Models
            2. Tuned Models
            3. Qualified Models
        """
        sns.set_theme(style="ticks", palette="Set3", font_scale=1.15)  # twilight_shifted_r
        # sns.set(style="ticks", font_scale=1.15)  # white, dark, whitegrid, darkgrid, ticks
        default_data = pd.DataFrame({
            "Experiment": [e for e in range(self.n_experiment)],
            "Acc_scores": self.scores["initial_scoring"]["acc_scores"],
            "Balan_acc_scores": self.scores["initial_scoring"]["balan_acc_scores"],
            "Models_type": "Default"
        })

        feature_selected_data = pd.DataFrame({
            "Experiment": [e for e in range(self.n_experiment)],
            "Acc_scores": self.scores["feature_selected_scoring"]["acc_scores"],
            "Balan_acc_scores": self.scores["feature_selected_scoring"]["balan_acc_scores"],
            "Models_type": "Default+FeatureSelected"
        })

        tuned_data = pd.DataFrame({
            "Experiment": [e for e in range(self.n_experiment)],
            "Acc_scores": self.scores["finalized_scoring"]["acc_scores"],
            "Balan_acc_scores": self.scores["finalized_scoring"]["balan_acc_scores"],
            "Models_type": "Tuned+FeatureSelected"
        })

        # Displaying `default_data`, `feature_selected_data` or `tuned_data`.
        display(tuned_data)

        print(tuned_data.describe())

        balanced_acc_threshold = tuned_data["Balan_acc_scores"].quantile(0.25)
        log.info(f"balanced_acc_threshold: {balanced_acc_threshold}")

        bad_models_ix = list(tuned_data[tuned_data["Balan_acc_scores"] < balanced_acc_threshold]["Experiment"])
        log.info(f"bad_models_ix: {bad_models_ix}")

        # remove bad_models_ix and add the another box plot
        qualified_models_ix = [e for e in range(self.n_experiment) if e not in bad_models_ix]
        log.info(f"qualified_models_ix: {qualified_models_ix}")
        tuned_qualified_data = tuned_data.iloc[qualified_models_ix].copy()
        tuned_qualified_data["Models_type"] = "Tuned+FeatureSelected+Qualified"

        # Assign qualified models
        self.qualified_models = [self.tuned_models[i] for i in qualified_models_ix]

        df_melted = pd.concat(
            [default_data, feature_selected_data, tuned_data, tuned_qualified_data]
        ).melt(
            id_vars=['Experiment', 'Models_type'],
            value_vars=['Acc_scores', 'Balan_acc_scores'],
            var_name='METRIC',
            value_name='SCORES'
        )

        # display(df_melted)

        comparison_data = pd.concat(
            [default_data, feature_selected_data, tuned_data, tuned_qualified_data]
        ).groupby('Models_type').mean().T.drop('Experiment')
        self.comparison_data = comparison_data

        sns.catplot(x='METRIC', y='SCORES', hue='Models_type', data=df_melted, kind=kind)
        # plt.axhline(y=0.7424242424242424, color='b', linestyle='--')
        # plt.axhline(y=0.6920995670995671, color='b', linestyle='--')

        plt.axhline(y=balanced_acc_threshold, color='r', linestyle='--')

        df_concated = pd.concat(
            [default_data, feature_selected_data, tuned_data, tuned_qualified_data]
        )
        log_simple.debug("{}\n".format(
            df_concated['Models_type'].value_counts().to_frame(name="Number of Model"))
        )


def cross_val_score_feature_comparison(X, y, scoring, n_repeats, n_jobs, clf=None):
    # In calculation of scores, cross-validation is repeated n times, which yields a total of 10*n folds.
    # E.g. if n=10, it means cross-validation is repeated 10 times with a total of 100 folds.
    if clf is None:
        clf = get_default_classifier(random_state=42)
    # Std of each scores can be added: scores.mean() → scores.std()
    return (round(cross_val_score(clf, X, y,
                                  cv=RepeatedStratifiedKFold(n_splits=10, n_repeats=n_repeats),
                                  scoring=scoring, n_jobs=n_jobs).mean(), 4))


def evaluate_metric(X_benchmark_feature_names_dataframes: dict, y, metric, n_repeats, n_jobs, verbose):
    scores_comparison = []
    for X_item_name, X_item in X_benchmark_feature_names_dataframes.items():
        scores = cross_val_score_feature_comparison(X_item, y, metric, n_repeats=n_repeats, n_jobs=n_jobs)
        scores_comparison.append(scores)
        if verbose:
            print("{: <28}: {}".format(X_item_name, scores))

    return scores_comparison


def evaluate_metrics(X_benchmark_feature_names_dataframes, y, n_repeats, n_jobs, verbose, eval_metrics_dict=None):
    if eval_metrics_dict is None:
        eval_metrics_dict = EVAL_METRICS_DICT

    scoring_metrics = {}

    for metric_name, metric in eval_metrics_dict.items():
        if verbose:
            print(F"\nEVALUATION METRIC: {metric_name}")
            print("------------------------------------")
        scores_comparison = evaluate_metric(X_benchmark_feature_names_dataframes, y, metric=metric,
                                            n_repeats=n_repeats,
                                            n_jobs=n_jobs, verbose=verbose)
        scoring_metrics[metric_name] = scores_comparison
        if verbose:
            print("====================================")

    scoring_metrics_table = pd.DataFrame(scoring_metrics, index=[feature_names for
                                                                 feature_names in X_benchmark_feature_names_dataframes])

    return scoring_metrics_table
