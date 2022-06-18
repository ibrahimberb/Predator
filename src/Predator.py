from typing import List, Tuple
from pathlib import Path

from tqdm.notebook import tqdm

from helpers.helpers_predator.data_sampling import prepare_data_spsm
from helpers.helpers_predator.paths import Paths
from helpers.helpers_predator.data_materials import DataMaterials
from helpers.helpers_predator.feature_selection import ShapFeatureSelector
from helpers.helpers_predator.evaluation import EvaluationMetrics, EvaluationValid
from helpers.helpers_predator.fine_tuning import FineTuner
from helpers.helpers_predator.predictions import predictions_object

from helpers.helpers_predator.common import (
    export_prediction_data,
    get_current_date_time,
)

from helpers.helpers_predator.models import (
    DefaultModels,
    TunedModels,
    QualifiedModels,
    FinalizedModels,
    EnsembleVotingClassifier
)

from helpers.mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)

TCGA_CODE = str


class Predator:
    def __init__(
            self,
            project_common_file_dir: Path,
            mutations_path: Path,
            initial_columns_path: Path,
            n_experiment: int,
            eliminate_models=False,
            random_seeds=None,
    ):

        log.debug('Initializing Predator ..')
        self.n_experiment = n_experiment
        self.n_models = self.n_experiment

        if random_seeds is None:
            self.random_seeds = list(range(1, self.n_experiment + 1))
        else:
            self.random_seeds = random_seeds

        self.paths = Paths(
            project_common_file_dir,
            mutations_path,
            initial_columns_path,
        )

        self.tcga_cohorts = []

        self.data_materials = DataMaterials(n_experiment, self.random_seeds)

        self.data_materials.initialize_train_dataset(
            project_common_file_dir, initial_columns_path, mutations_path
        )

        self.default_models = DefaultModels(self.n_experiment)
        self.tuned_models = None
        self.finalized_models = None
        self.qualified_models = None

        self.eval_valid = EvaluationValid(
            self.n_experiment,
            self.data_materials,
            self.default_models,
            self.tuned_models,
        )

        self.determined_feature_set = None
        self.determined_features = None

        self.shap_feature_selector = None
        self.eval_metrics = None
        self.fine_tuner = None

        self.ensemble_voting_classifier = None
        self.predictions = None

        self.eliminate_models = eliminate_models

        self.config = {}

        # variant = spsm
        # self.variant = variant

    def sample_spsm(self):
        log.debug("sampling ..")
        sampled_train_data_list = []
        for i in tqdm(range(self.n_experiment)):
            sampled_train_data = prepare_data_spsm(
                self.data_materials["train_data_processed"], random_seed=self.random_seeds[i]
            )
            sampled_train_data_list.append(sampled_train_data)

        self.data_materials["sampled_train_data_list"] = sampled_train_data_list

    def run_evaluate_valid(
        self,
        models_type="default",
        show_confusion_matrix=False,
    ):

        if models_type == 'tuned':
            self.eval_valid.tuned_models = self.tuned_models

        self.eval_valid.evaluate(
            models_type=models_type,
            show_confusion_matrix=show_confusion_matrix,
            determined_features=self.determined_features
        )

    def init_shap_feature_selector(self, shap_top_ns):
        self.shap_feature_selector = ShapFeatureSelector(self.n_experiment, shap_top_ns)
        self.shap_feature_selector.load_shap_values(self.data_materials)  # fixme -> select_features
        self.shap_feature_selector.get_selected_features(self.data_materials)

    def aggregate_selected_features(self, method):
        self.shap_feature_selector.shap_aggregate_selected_features(method)
        for shap_top_n, aggregated_feature in self.shap_feature_selector.n_features_to_aggregated_features.items():
            self.data_materials.initialize_feature_selected_data_materials(shap_top_n)
            self.data_materials.append_feature_selected_data_materials(shap_top_n, aggregated_feature)

    def initialize_evaluation_metrics(self):
        self.eval_metrics = EvaluationMetrics(self.n_experiment, self.data_materials, self.shap_feature_selector)

    def set_determined_feature_set(self, feature_set: str):
        _, top_n = feature_set.split('_')
        determined_features = (
            self
            .shap_feature_selector
            .n_features_to_aggregated_features[int(top_n)]
        )
        log.debug(f"Setting determined feature set to `{feature_set}`.")
        self.determined_feature_set = feature_set
        log.debug(f"Setting determined features to \n{determined_features}.")
        self.determined_features = determined_features

    # def init_fine_tuner(self, n_iter, n_repeats_cv, n_jobs, verbose):
    #     # Fine tuning will be applied to training set, not all training data.
    #     Xs_determined = f"Xs_train_{self.determined_feature_set}"
    #
    #     self.fine_tuner = FineTuner(
    #         n_iter=n_iter, n_repeats_cv=n_repeats_cv, n_jobs=n_jobs, verbose=verbose,
    #         Xs_determined=Xs_determined, n_experiment=self.n_experiment,
    #         random_seeds=self.random_seeds, data_materials=self.data_materials
    #     )

    def run_hyperparameter_search(
            self,
            n_iter,
            n_repeats_cv,
            n_jobs,
            verbose,
            search_type,
            param_grid_level,
    ):

        # Fine tuning will be applied to training set, not all training data.
        Xs_determined = f"Xs_train_{self.determined_feature_set}"

        self.fine_tuner = FineTuner(
            n_iter=n_iter,
            n_repeats_cv=n_repeats_cv,
            n_jobs=n_jobs, verbose=verbose,
            Xs_determined=Xs_determined,
            n_experiment=self.n_experiment,
            random_seeds=self.random_seeds,
            data_materials=self.data_materials,
            param_grid_level=param_grid_level
        )

        hyperparam_config = {
            "n_iter": n_iter,
            "n_repeats_cv": n_repeats_cv,
            "search_type": search_type,
            "param_grid": self.fine_tuner.param_grid
        }

        self.config['hyperparameters'] = hyperparam_config

        self.fine_tuner.run_search(search_type)
        self.tuned_models = TunedModels(self.fine_tuner.best_estimators)

    def compare_tuned_models(self, kind='box'):
        self.eval_valid.compare_models(kind)
        self.qualified_models = QualifiedModels(self.eval_valid.qualified_models)

    def fit_finalized_models(self):
        log.debug("Fitting finalized models with all training data ..")
        if self.eliminate_models:
            log.info(f"Model elimination: {self.eliminate_models}")
            log.info(f"Using {len(self.qualified_models)} qualified models as finalized models.")
            self.finalized_models = FinalizedModels(self.qualified_models)
            self.n_models = len(self.qualified_models)

        else:
            log.info(f"Model elimination: {self.eliminate_models}")
            log.info(f"Using {len(self.tuned_models)} tuned models as finalized models.")
            self.finalized_models = FinalizedModels(self.tuned_models)

        # Fill the config main
        config_main = {
            "eliminate_models": self.eliminate_models,
            "n_experiment": self.n_experiment,
            "n_models": self.n_models,
            "preparation_completed": get_current_date_time()
        }
        self.config["main"] = config_main

        # Fill each model
        self.finalized_models.fit_all(self.data_materials, self.determined_feature_set)

    def initialize_target_data_materials(
            self,
            tcga_code_path_pairs: List[Tuple[TCGA_CODE, Path]],
            initial_columns_path=None
    ):

        self.tcga_cohorts = []
        for tcga, _ in tcga_code_path_pairs:
            self.tcga_cohorts.append(tcga)

        if initial_columns_path is None:
            initial_columns_path = self.paths.initial_columns_path

        self.data_materials.initialize_target_datasets(
            self.paths.project_common_file_dir,
            initial_columns_path,
            tcga_code_path_pairs
        )

        self.data_materials.initialize_target_data_materials(
            self.determined_features, tcga_code_path_pairs
        )

    def predict(self, voting="hard"):
        log.debug("Predicting on cancer datasets ..")
        self.ensemble_voting_classifier = EnsembleVotingClassifier(
            self.finalized_models, voting=voting
        )

        self.predictions = predictions_object(voting, self.n_models)

        for tcga in self.tcga_cohorts:
            log.debug(f"Predicting on {tcga} cohort ..")
            tcga_predictions = self.ensemble_voting_classifier.predict(
                self.data_materials[f"Xs_{tcga}"]
            )
            self.predictions.add_predictions(tcga, tcga_predictions)

    def predictions_post_process(self):
        for tcga in self.tcga_cohorts:
            self.predictions.post_process_predictions(
                tcga,
                self.data_materials[f"{tcga}"]
            )

    def prepare_ensemble_prediction_data(self):
        for tcga in self.tcga_cohorts:
            self.predictions.prepare_ensemble_prediction_data(
                tcga, self.data_materials[f"{tcga}"]
            )

    def export_prediction(
        self,
        tcga,
        data,
        file_name,
        folder_path,
        voting,
        overwrite=False,
        file_extension='csv'
    ):

        # Fill the config main
        config_prediction = {
            # TODO: add info about tcga data: (its path/name and dimensions)
            "tcga_cohorts": self.tcga_cohorts,
        }
        self.config["prediction"] = config_prediction

        export_prediction_data(
            tcga=tcga,
            data=data,
            file_name=file_name,
            folder_path=folder_path,
            config=self.config,
            overwrite=overwrite,
            voting=voting,
            file_extension=file_extension
        )
