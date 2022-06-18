from typing import List

import seaborn as sns

import pandas as pd
from pandas import DataFrame
from ..mylogger import get_handler
from .predictions_utils import (
    get_predictive_columns_removed_data,
    max_votes,
    convert_primary_isomer,
    get_triplet_columns,
    remove_triplet_columns,
    drop_invalid_predicted_entries,
    drop_invalid_predicted_probs_entries,
    add_votes,
    add_voted_probs,
    take_avg,
    take_median,
)

import logging
import matplotlib.pyplot as plt

from ..labels import ClassLabels

from tqdm.notebook import tqdm

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


def predictions_object(voting, n_models):
    if voting == 'hard':
        return PredictionsHard(n_models)
    elif voting == 'soft':
        return PredictionsSoft(n_models)


class Predictions(dict):
    def __init__(self, n_models):
        super().__init__()
        self.n_models = n_models
        self.value_counts = {}
        self.predictions_distributions_per_model = {}

    def init_value_counts(self, tcga):
        log.debug("Initializing value counts ..")
        value_counts = []
        for model in range(self.n_models):
            pred = self[tcga][model]
            value_counts.append((len(pred[pred == ClassLabels.DISRUPTING]), len(pred[pred == ClassLabels.NONDISRUPTING])))
        self.value_counts[tcga] = value_counts

    def prepare_predictions_distribution_data(self, tcga):
        predictions_distributions = pd.DataFrame(
            self.value_counts[tcga],
            index=[f"TRIAL_{model}" for model in range(1, self.n_models + 1)],
        )
        predictions_distributions.rename(
            {
                ClassLabels.DISRUPTING : 'Disrupting', 
                ClassLabels.NONDISRUPTING: 'NoEffect+Increasing'
            },
            axis="columns",
            inplace=True
        )

        predictions_distributions["Total_entry"] = (
                predictions_distributions["Disrupting"]
                + predictions_distributions["NoEffect+Increasing"]
        )

        self.predictions_distributions_per_model[tcga] = predictions_distributions

    def plot_distributions(self, tcga):
        sns.set(style="white", font_scale=1.15)
        predictions_distributions_data = self.predictions_distributions_per_model[tcga].copy()
        predictions_distributions_data = (
            predictions_distributions_data
            .rename_axis("MODEL")
            .reset_index()
        )
        predictions_distributions_data.plot(
            x="MODEL",
            y=["Label_X", "Label_Y"],    # ["Disrupting", "NoEffect+Increasing"],
            kind="bar",
            figsize=(15, 5),
        )
        plt.grid(zorder=0, axis="y")
        plt.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
        plt.title("Distribution of Class Predictions per Model")
        plt.xlabel("Models")
        plt.ylabel("Counts")
        plt.show()

    def plot_predictions_distributions(self, tcga):
        self.init_value_counts(tcga)
        self.prepare_predictions_distribution_data(tcga)
        self.plot_distributions(tcga)

    def plot_distribution_valid_vs_invalid(self, tcga):
        valid_vs_invalid_counts_data = pd.DataFrame({
            "Num_Valid": [len(data) for data in self[f"{tcga}_predicted_valid_datasets"]],
            "Num_Invalid": [len(data) for data in self[f"{tcga}_predicted_invalid_datasets"]]
        })

        valid_vs_invalid_counts_data.index = [f'TRIAL_{i}' for i in range(1, self.n_models + 1)]
        valid_vs_invalid_counts_data.index.name = 'MODEL'
        valid_vs_invalid_counts_data.reset_index().plot(
            x="MODEL", y=["Num_Valid", "Num_Invalid"], kind="bar",
            figsize=(15, 5), color=['#9bb7bf', '#ed3e4f']
        )
        plt.grid(zorder=0, axis='y')
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        plt.title(f"{tcga.upper()} Prediction Distribution of Number of Valid Entries vs Invalid (Dropped) Entries")
        plt.xlabel("Models")
        plt.ylabel("Number of Entries")
        plt.show()

    def plot_num_finalized_predictions(self, tcga):
        log.debug("Plotting number of finalized predictions per model.\n"
                  "Note that following plot shows the number of (protein, mutation, interactor) "
                  "triplets which had valid prediction.")
        finalized_prediction_counts_data = pd.DataFrame({
            "Num_Entries": [len(data) for data in self[f"{tcga}_finalized_prediction_dataframes"]],
        })

        finalized_prediction_counts_data.index = [f'TRIAL_{i}' for i in range(1, self.n_models + 1)]
        finalized_prediction_counts_data.index.name = 'MODEL'
        finalized_prediction_counts_data.reset_index().plot(
            x="MODEL", y=["Num_Entries"], kind="bar",
            figsize=(15, 5), color=['#24BFA5']
        )
        plt.grid(zorder=0, axis='y')
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        plt.title(f"{tcga.upper()} Number of Unique Entries per Model in Finalized Prediction Dataframes")
        plt.xlabel("Models")
        plt.ylabel("Number of Entries")
        plt.show()

    def get_predictive_columns_removed_datasets(
            self,
            prediction_datasets: List[DataFrame]
    ) -> List[DataFrame]:
        """
        Remove predictive columns (i.e. feature columns) from each datasets, leaving them
        with predictions and triplets.

        Parameters
        ----------
            prediction_datasets : <List[DataFrame]>
                A list of dataframe containing cancer data with predictions.

        Returns
        -------
            features_removed_datasets : <List[DataFrame]>
                List of dataframes, each containing following columns:
                ["Prediction", "UniProt_ID", "Mutation", "Interactor_UniProt_ID"]
        """

        features_removed_datasets = []
        for m in range(self.n_models):
            features_removed_data = get_predictive_columns_removed_data(prediction_datasets[m])
            features_removed_datasets.append(features_removed_data)

        return features_removed_datasets

    def prepare_finalized_prediction_datasets(self, tcga, prediction_datasets):
        tcga_finalized_prediction_dataframes = self.get_predictive_columns_removed_datasets(
            prediction_datasets
        )
        self.drop_duplicated_entries_datasets(tcga_finalized_prediction_dataframes)
        self[f"{tcga}_finalized_prediction_dataframes"] = tcga_finalized_prediction_dataframes

    def drop_duplicated_entries_datasets(self, datasets: List[DataFrame]) -> None:
        """
        Drops the duplicated entries in each dataframes.
        Resets the index after dropping entries.

        Parameters
        ----------
            datasets : <List[DataFrame]>

        Returns
        -------
            None, modifies the dataframes inplace.
        """
        for m in range(self.n_models):
            datasets[m].drop_duplicates(keep="first", inplace=True)
            datasets[m].reset_index(drop=True, inplace=True)

    def plot_ensemble_prediction_distribution(self, tcga):
        voted_predictions = self[f"{tcga}_ensemble_prediction_data"]['VOTED_PREDICTION']
        voted_predictions.value_counts().plot(kind="bar")
        plt.title(f'Distribution of predictions for {tcga}')
        plt.xlabel('Predictions')
        plt.ylabel('Count')
        plt.grid()


class PredictionsHard(Predictions):
    def __init__(self, n_models):
        log.debug(f'Initializing: {__class__.__name__}')
        super().__init__(n_models=n_models)

    def add_predictions(self, tcga, tcga_predictions):
        log.debug(f"{__class__.__name__}")
        log.debug(f"Predicting on {tcga} cohort ..")
        log.debug(f"Adding key `{tcga}` to self.predictions")
        self[tcga] = tcga_predictions

    def prepare_ensemble_prediction_data(self, tcga, tcga_data):
        log.debug(f"{__class__.__name__}")
        log.debug(f"Preparing ensemble prediction data for {tcga} ..")
        tcga_ensemble_prediction_data = get_triplet_columns(tcga_data)
        tcga_ensemble_prediction_data = convert_primary_isomer(
            "Interactor_UniProt_ID", tcga_ensemble_prediction_data
        )
        tcga_ensemble_prediction_data.drop_duplicates(keep="first", inplace=True)
        tcga_ensemble_prediction_data.reset_index(drop=True, inplace=True)
        # brca_ensemble_prediction_data.apply(lambda voting: voting()) ## todo vectorication using 3 column as param.
        tcga_ensemble_prediction_data = add_votes(
            tcga_ensemble_prediction_data, self[f"{tcga}_finalized_prediction_dataframes"]
        )
        tcga_ensemble_prediction_data['VOTED_PREDICTION'] = tcga_ensemble_prediction_data.apply(
            max_votes, axis=1
        )

        self[f"{tcga}_ensemble_prediction_data"] = tcga_ensemble_prediction_data
        log.debug(f"Ensemble prediction data for {tcga} is prepared.")

        tcga_prediction_results = tcga_ensemble_prediction_data.drop(
            [f"Num_preds_{ClassLabels.DISRUPTING}", f"Num_preds_{ClassLabels.NONDISRUPTING}", "Num_preds_NO_VOTE"], axis='columns'
        )

        # Rename "VOTED_PREDICTION" column to "Prediction"
        tcga_prediction_results = tcga_prediction_results.rename(
            columns={'VOTED_PREDICTION': 'Prediction'}
        )
        self[f"{tcga}_prediction_results"] = tcga_prediction_results
        log.debug(f"Resulting prediction data is available for {tcga}.\n"
                  f"Accessible from predictions.['{tcga}_prediction_results']")

        tcga_prediction_results_no_votes_dropped = (
            tcga_prediction_results[tcga_prediction_results['Prediction'].isin([ClassLabels.DISRUPTING, ClassLabels.NONDISRUPTING])]
        )
        self[f"{tcga}_prediction_results_no_votes_dropped"] = tcga_prediction_results_no_votes_dropped
        log.debug(f"Resulting prediction data (no_votes dropped) is available for {tcga}.\n"
                  f"Accessible from predictions.['{tcga}_prediction_results_no_votes_dropped']")

    def merge_predictions_cancer_datasets(self, tcga, tcga_data):
        log.debug(f"Merging predictions with {tcga} cancer dataset ..")
        tcga_predicted_datasets = []
        for m in range(self.n_models):
            tcga_predicted_data = tcga_data.copy(deep=True)
            # Insert prediction array values into data as first (0th) column.
            tcga_predicted_data.insert(0, 'Prediction', self[f"{tcga}"][m])
            tcga_predicted_datasets.append(tcga_predicted_data)

        self[f"{tcga}_predicted_datasets"] = tcga_predicted_datasets

    def post_process_predictions(self, tcga, tcga_datasets):
        log.debug(f'{__class__.__name__}')
        log.debug(f"Post processing predictions for cohort {tcga} ..")

        self.merge_predictions_cancer_datasets(
            tcga, tcga_datasets
        )

        log.debug(f"Handling valid and invalid entries ..")
        tcga_predicted_valid_datasets = []
        tcga_predicted_invalid_datasets = []
        for m in tqdm(range(self.n_models)):
            isomer_converted_data = convert_primary_isomer(
                column_name="Interactor_UniProt_ID",
                data=self[f"{tcga}_predicted_datasets"][m]
            )
            predicted_valid_data, predicted_invalid_data = drop_invalid_predicted_entries(
                isomer_converted_data
            )

            tcga_predicted_valid_datasets.append(predicted_valid_data)
            tcga_predicted_invalid_datasets.append(predicted_invalid_data)

        self[f"{tcga}_predicted_valid_datasets"] = tcga_predicted_valid_datasets
        self[f"{tcga}_predicted_invalid_datasets"] = tcga_predicted_invalid_datasets
        log.debug(f"Preparing finalized prediction datasets for {tcga} ..")
        self.prepare_finalized_prediction_datasets(tcga, tcga_predicted_valid_datasets)
        log.debug(f"Post processing completed for {tcga}.")


class PredictionsSoft(Predictions):
    _note = "used_class_1_prob"

    def __init__(self, n_models):
        log.debug(f'Initializing: {__class__.__name__}')
        super().__init__(n_models=n_models)

    def add_predictions(self, tcga, tcga_predictions_probabilities):
        log.debug(f'{__class__.__name__}')
        log.debug(f"Predicting probabilities on {tcga} cohort ..")
        log.debug(f"Adding key `{tcga}_prob` to self.predictions")
        self[f"{tcga}_prob"] = tcga_predictions_probabilities

    def merge_predictions_cancer_datasets(self, tcga, tcga_data):
        log.debug(f'{__class__.__name__}')
        log.debug(f"Merging predictions with {tcga} cancer dataset ..")
        tcga_predicted_probs_datasets = []
        for m in range(self.n_models):
            tcga_predicted_probs_data = tcga_data.copy(deep=True)
            # Insert predictions probability of class 1 (whatever that might be) array values into data as first (0th) column.
            tcga_predicted_probs_data.insert(0, 'Prediction', self[f"{tcga}_prob"][m][:, 1])
            tcga_predicted_probs_datasets.append(tcga_predicted_probs_data)

        self[f"{tcga}_predicted_probs_datasets"] = tcga_predicted_probs_datasets

    def post_process_predictions(self, tcga, tcga_datasets):
        log.debug(f'{__class__.__name__}')
        log.debug(f"Post processing predictions for cohort {tcga} ..")

        self.merge_predictions_cancer_datasets(
            tcga, tcga_datasets
        )

        log.debug(f"Handling valid and invalid entries ..")
        tcga_predicted_valid_datasets = []
        tcga_predicted_invalid_datasets = []
        for m in tqdm(range(self.n_models)):
            isomer_converted_data = convert_primary_isomer(
                column_name="Interactor_UniProt_ID",
                data=self[f"{tcga}_predicted_probs_datasets"][m]
            )
            predicted_valid_data, predicted_invalid_data = drop_invalid_predicted_probs_entries(
                isomer_converted_data
            )

            tcga_predicted_valid_datasets.append(predicted_valid_data)
            tcga_predicted_invalid_datasets.append(predicted_invalid_data)

        self[f"{tcga}_predicted_valid_datasets"] = tcga_predicted_valid_datasets
        self[f"{tcga}_predicted_invalid_datasets"] = tcga_predicted_invalid_datasets
        log.debug(f"Preparing finalized prediction datasets for {tcga} ..")
        self.prepare_finalized_prediction_datasets(tcga, tcga_predicted_valid_datasets)
        log.debug(f"Post processing completed for {tcga}.")

    def prepare_ensemble_prediction_data(self, tcga, tcga_data, take="median"):
        # FIXME: REFORMATTING ...
        log.info(f'{__class__.__name__}')
        log.info(f"Preparing ensemble prediction data for {tcga} taking {take} ..")
        tcga_ensemble_prediction_data = get_triplet_columns(tcga_data)
        tcga_ensemble_prediction_data = convert_primary_isomer(
            "Interactor_UniProt_ID", tcga_ensemble_prediction_data
        )

        tcga_ensemble_prediction_data.drop_duplicates(keep="first", inplace=True)
        tcga_ensemble_prediction_data.reset_index(drop=True, inplace=True)

        tcga_ensemble_prediction_data = add_voted_probs(
            tcga_ensemble_prediction_data, self[f"{tcga}_finalized_prediction_dataframes"]
        )

        tcga_predictions_prob_data = remove_triplet_columns(tcga_ensemble_prediction_data)

        if take == "average":
            tcga_predictions_prob_data['PROB_1s_AVG'] = tcga_predictions_prob_data.apply(
                take_avg, axis=1
            )
        elif take == "median":
            tcga_predictions_prob_data['PROB_1s_AVG'] = tcga_predictions_prob_data.apply(
                take_median, axis=1
            )
        else:
            raise ValueError(f"Invalid value for parameter `take`: {take}")

        # If probability of being class 1 is greater than or equal to 0.50,
        # it is assigned as class 1.
        # tcga_predictions_prob_data['VOTED_PREDICTION'] = (
        #     (tcga_predictions_prob_data['PROB_1s_AVG'] > 0.50).astype('int')
        # )

        # Unless PROB_1s_AVG is NO_VOTE, if probability of being class 1 is greater than 0.50,
        # it is assigned as class 1.
        tcga_predictions_prob_data['VOTED_PREDICTION'] = tcga_predictions_prob_data['PROB_1s_AVG'].apply(
            lambda x: 'NO_VOTE' if x == 'NO_VOTE' else int(float(x) > 0.50)
        )

        self[f"{tcga}_predictions_prob_data"] = tcga_predictions_prob_data
        log.debug(f"Prediction probabilities data for {tcga} is prepared.\n"
                  f"Accessible from `{tcga}_predictions_prob_data`.")

        tcga_ensemble_prediction_data['VOTED_PREDICTION'] = tcga_predictions_prob_data['VOTED_PREDICTION']

        self[f"{tcga}_ensemble_prediction_data"] = tcga_ensemble_prediction_data
        log.debug(f"Ensemble prediction data for {tcga} is prepared.\n"
                  f"Accessible from `{tcga}_ensemble_prediction_data`.")

        tcga_prediction_results = get_triplet_columns(tcga_ensemble_prediction_data)
        tcga_prediction_results['Prediction'] = tcga_ensemble_prediction_data['VOTED_PREDICTION']
        self[f"{tcga}_prediction_results"] = tcga_prediction_results
        log.debug(f"Resulting prediction data is available for {tcga}.\n"
                  f"Accessible from predictions.['{tcga}_prediction_results']")

        tcga_prediction_results_no_votes_dropped = (
            tcga_prediction_results[tcga_prediction_results['Prediction'].isin([0, 1])]
        )
        self[f"{tcga}_prediction_results_no_votes_dropped"] = tcga_prediction_results_no_votes_dropped
        log.debug(f"Resulting prediction data (no_votes dropped) is available for {tcga}.\n"
                  f"Accessible from predictions.['{tcga}_prediction_results_no_votes_dropped']")
