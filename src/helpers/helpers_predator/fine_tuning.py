from typing import List

from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import RepeatedStratifiedKFold

import pandas as pd
import numpy as np

from .machine_learning_utils import get_default_classifier

from tqdm.notebook import tqdm

from ..mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)

# PARAM_GRID_SIMPLE = {
#     'bootstrap': [True, False],
#     'max_depth': [2, 5, 10],
#     'max_features': ['auto', 'sqrt'],
#     'n_estimators': [50, 100, 200, 400],
#     'min_samples_split': [2, 5, 10]
# }


class FineTuner:

    PARAM_GRID_LEVEL_0 = {
        'max_depth': [2, 5, 10],
        'n_estimators': [10, 25, 50, 75, 100, 200, 400],
        'min_samples_split': [2, 5],
        'max_features': ['auto', 'sqrt', None],
        'class_weight': ['balanced', None]
    }

    PARAM_GRID_LEVEL_1 = {
        'bootstrap': [True, False],
        'max_depth': list(np.arange(2, 15)) + [None],
        'min_samples_leaf': [1, 2, 4],
        'n_estimators': [5] + list(np.arange(10, 110, 10)) + list(np.arange(120, 620, 20)),
        'min_samples_split': [2, 5],
        'max_features': ['auto', 'sqrt', None],
        'class_weight': ['balanced', None]
    }

    def __init__(
            self,
            n_iter: int,
            n_repeats_cv: int,
            n_jobs: int,
            verbose: int,
            Xs_determined,
            n_experiment: int,
            random_seeds: List[int],
            data_materials,
            param_grid_level=0
    ):
        self.n_iter = n_iter
        self.n_repeats_cv = n_repeats_cv
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.n_experiment = n_experiment
        self.random_seeds = random_seeds
        self.data_materials = data_materials
        self.Xs_determined = Xs_determined

        self.randomized_search_objects = None
        self.classifiers_attributes_data = None
        self.best_estimators = None

        self.param_grid = self.get_param_grid(param_grid_level)

    def get_param_grid(self, level):
        if level == 0:
            return self.PARAM_GRID_LEVEL_0

        elif level == 1:
            return self.PARAM_GRID_LEVEL_1

        else:
            raise ValueError(f"Level {level} is invalid argument.")

    # def run_randomized_search(self) -> None:
    #     log.debug("Running randomized search for each experiment ..")
    #     randomized_search_objects = []
    #     for exp in tqdm(range(self.n_experiment)):
    #         random_seed = self.random_seeds[exp]
    #         X_fine_tuning = self.data_materials[self.Xs_determined][exp]
    #         y_fine_tuning = self.data_materials["ys_train"][exp]
    #         randomized_search = self.get_randomized_search(random_seed)
    #         randomized_search.fit(X_fine_tuning, y_fine_tuning)
    #         randomized_search_objects.append(randomized_search)
    #     self.randomized_search_objects = randomized_search_objects
    #     self.save_randomized_search_info()
    #     self.save_best_estimators()

    def run_search(self, search_type="randomized") -> None:
        log.debug(f"Running {search_type} search for each experiment ..")
        log.debug(f"PARAM_GRID: {self.param_grid}")
        search_objects = []
        # searcher = self.get_randomized_search if search_type == "randomized" else self.get_grid_search
        for exp in tqdm(range(self.n_experiment)):
            random_seed = self.random_seeds[exp]
            X_fine_tuning = self.data_materials[self.Xs_determined][exp]
            y_fine_tuning = self.data_materials["ys_train"][exp]
            search_obj = self.get_searcher(search_type, random_seed)
            search_obj.fit(X_fine_tuning, y_fine_tuning)
            search_objects.append(search_obj)
        self.randomized_search_objects = search_objects
        self.save_randomized_search_info()
        self.save_best_estimators()

    def get_searcher(self, search_type, random_seed):
        if search_type == "randomized":
            return self.get_randomized_search(random_seed)
        elif search_type == "grid":
            return self.get_grid_search(random_seed)
        else:
            raise ValueError(f"{search_type} is invalid argument.")

    def get_randomized_search(self, random_seed) -> RandomizedSearchCV:

        if self.n_iter is None:
            raise ValueError("Using RANDOMIZED SEARCH. Parameter `n_iter` cannot be None.")

        clf = get_default_classifier(random_state=random_seed)

        randomized_search = RandomizedSearchCV(
            clf, self.param_grid,
            n_iter=self.n_iter,
            random_state=random_seed,
            cv=RepeatedStratifiedKFold(
                n_splits=10, n_repeats=self.n_repeats_cv, random_state=random_seed
            ),
            scoring='balanced_accuracy',
            return_train_score=True,
            n_jobs=self.n_jobs,
            verbose=self.verbose
        )

        return randomized_search

    def get_grid_search(self, random_seed) -> GridSearchCV:

        if self.n_iter is not None:
            log.warning("Using GRID SEARCH. Parameter `n_iter` will be ignored.")

        clf = get_default_classifier(random_state=random_seed)

        grid_search = GridSearchCV(
            clf,
            self.param_grid,
            cv=RepeatedStratifiedKFold(
                n_splits=10, n_repeats=self.n_repeats_cv, random_state=random_seed
            ),
            scoring='balanced_accuracy',
            return_train_score=True,
            n_jobs=self.n_jobs,
            verbose=self.verbose
        )

        return grid_search

    def save_randomized_search_info(self):
        experiment_repeat_to_randomized_search_info = {}

        for exp in range(self.n_experiment):
            experiment_repeat_to_randomized_search_info[F'EXP_{exp + 1}'] = [
                self.randomized_search_objects[exp].best_params_,
                self.randomized_search_objects[exp].best_estimator_,
                self.randomized_search_objects[exp].best_score_]

        classifiers_attributes_data = pd.DataFrame(
            experiment_repeat_to_randomized_search_info,
            index=['best_params_', 'best_estimator_', 'best_score_']
        ).T

        self.classifiers_attributes_data = classifiers_attributes_data

    def save_best_estimators(self):
        best_estimators = []
        for search_obj in self.randomized_search_objects:
            best_estimators.append(search_obj.best_estimator_)
        self.best_estimators = best_estimators
