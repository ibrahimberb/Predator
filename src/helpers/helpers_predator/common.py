import pickle

from pandas import DataFrame
import numpy as np
import os.path as op
from datetime import datetime
from pathlib import Path

from ..mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


def print_annotation(s):
    print(f"\n{s}\n{'-' * len(s)}")


def get_file_path(project_common_file_dir, filename):
    return op.join(Path(project_common_file_dir), Path(filename))


def export_data(
        folder_path: Path,
        data: DataFrame,
        file_name: str,
        file_extension='csv',
        overwrite=False
) -> None:
    """
    A helper function to export given data with specified name and extension.
    """

    log.debug(f"Exporting data {file_name} at location {folder_path} ..")
    file_name = op.join(folder_path, file_name)
    file_date = datetime.today().strftime('%Y-%m-%d')
    file_name = f'{file_name}_{file_date}.{file_extension}'

    # Ensure the file is not exists before creating to prevent overwriting.
    if op.isfile(file_name) and not overwrite:
        log.error(f"File {file_name} is already exist.\n"
                  "To overwrite existing file, use `overwrite=True`.")

    else:
        # Export
        data.to_csv(file_name, index=False)
        log.info(f'{file_name} is exported successfully.')


def get_random_id(path):
    import uuid

    while True:
        random_id = uuid.uuid4().hex[:8]
        id_path = op.join(path, random_id)
        is_valid = not op.isdir(id_path)

        if is_valid:
            log.debug(f"Folder with ID {random_id} is created.")
            return random_id


def export_serialized_predator(predator, folder_path="PredatorModels") -> None:
    """
    Serialize the `predator` object and export it as PKL file.

    Parameters
    ----------
        predator : <Predator>
            The predator object to be serialized and exported.

        folder_path : <PredatorModels>
            The folder in which predator models are located at.

    Returns
    -------
        None. Exports to file.
    """
    # TODO: Create a folder named `SerializedModel` (find better name) and
    #  export the config file, too.
    #  folder will be containing production ready pickle object and its metadata.

    current_date = datetime.today().strftime('%Y-%m-%d')

    random_id = get_random_id(folder_path)
    folder_name = op.join(f"PredatorModel_{current_date}", random_id)
    log.debug(f"Exporting Predator at location {folder_path} in folder {folder_name}..")

    folder_path = op.join(folder_path, folder_name)
    Path(f"{folder_path}").mkdir(parents=True, exist_ok=True)

    filename = op.join(folder_path, "predator.pkl")

    # Ensure the file is not exists before creating to prevent overwriting.
    if op.isfile(filename):
        log.error(f"File {filename} is already exist.\n"
                  "To overwrite existing file, use `overwrite=True`.")

    with open(filename, 'wb') as fid:
        pickle.dump(predator, fid)

    log.info(f"Predator object {filename} is exported.")

    # Export metadata
    config_path = op.join(folder_path, "config.txt")
    export_config(predator.config, config_path)


def load_predator(predator_obj_name):
    with open(predator_obj_name, 'rb') as inp:
        predator = pickle.load(inp)

    log.info(f"Predator object {predator_obj_name} is loaded successfully.")
    return predator


def export_prediction_data(
        tcga,
        data,
        file_name,
        folder_path,
        config,
        overwrite,
        voting,
        file_extension='csv'
) -> None:
    """
    A helper function to export given data with specified name and extension along with metadata.
    """
    current_date = datetime.today().strftime('%Y-%m-%d')
    current_time = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

    random_id = get_random_id(folder_path)
    folder_name = op.join(f"{tcga}_prediction_{current_date}", random_id)
    log.debug(f"Exporting data {file_name} at location {folder_path} in folder {folder_name}..")

    folder_path = op.join(folder_path, folder_name)
    file_name = op.join(folder_path, file_name)
    file_name = f"{file_name}_{voting}_{current_date}.{file_extension}"

    Path(f"{folder_path}").mkdir(parents=True, exist_ok=True)

    # Ensure the file is not exists before creating to prevent overwriting.
    if op.isfile(file_name) and not overwrite:
        log.error(f"File {file_name} is already exist.\n"
                  "To overwrite existing file, use `overwrite=True`.")

    else:
        # Export data
        data.to_csv(file_name, index=False)
        log.info(f'{file_name} is exported successfully.')

        # Export metadata
        config_path = op.join(folder_path, "config.txt")
        export_config(config, config_path)


def export_config(config, config_path):
    with open(config_path, "w") as file:
        for config_name, config_dict in config.items():
            file.write(f"{config_name.upper()}\n")
            for key, values in config_dict.items():
                file.write(
                    f"    {key.upper()}: {config[f'{config_name}'][f'{key}']}\n"
                )
    log.info(f"Config is exported.")


def get_current_date_time():
    current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    return current_time


def compare_predator_objects(p1, p2):
    """
    Compares given Predator objects.

    Parameters
    ----------
        p1 : <Predator>
            The first predator object.

        p2 : <Predator>
            The first predator object.

    Returns
    -------
        None if given Predators are equal.
    """

    assert p1.__dict__.keys() == p2.__dict__.keys()

    for k in p1.__dict__.keys():
        log.debug(f"Checking attribute {k} ..")

        if k in ["n_experiment", "n_models", "random_seeds", "tcga_cohorts", "eliminate_models", "config"]:
            assert getattr(p1, k) == getattr(p2, k)

        elif k == "paths":
            assert p1.paths.__dict__.keys() == p2.paths.__dict__.keys()
            for k_paths in p1.paths.__dict__.keys():
                assert getattr(p1.paths, k_paths) == getattr(p2.paths, k_paths)

        elif k == "data_materials":
            assert p1.data_materials.__dict__.keys() == p2.data_materials.__dict__.keys()
            assert p1.data_materials.keys() == p2.data_materials.keys()

            for k_data_materials in p1.data_materials.__dict__.keys():
                assert getattr(p1.data_materials, k_data_materials) == \
                       getattr(p2.data_materials, k_data_materials)

            for k_data_materials in p1.data_materials.keys():
                # print(f"datamaterial - {k_data_materials}")
                if isinstance(p1.data_materials[k_data_materials], list):
                    for k_ix, _ in enumerate(p1.data_materials[k_data_materials]):
                        assert p1.data_materials[k_data_materials][k_ix].equals(
                            p2.data_materials[k_data_materials][k_ix])

                elif isinstance(p1.data_materials[k_data_materials], DataFrame):
                    assert p1.data_materials[k_data_materials].equals(
                        p2.data_materials[k_data_materials])

        elif k in ["default_models", "tuned_models", "finalized_models", "qualified_models"]:
            assert len(p1.default_models) == len(p2.default_models)
            for ix, _ in enumerate(p1.default_models):
                assert np.array_equal(
                    p1.default_models[ix].classes_, p2.default_models[ix].classes_
                )

                assert str(p1.default_models[ix]) == str(p2.default_models[ix])

        elif k == "eval_valid":
            assert p1.eval_valid.__dict__.keys() == p2.eval_valid.__dict__.keys()
            assert p1.eval_valid.scores == p2.eval_valid.scores
            p1.eval_valid.comparison_data.equals(p2.eval_valid.comparison_data)

        elif k in ["determined_feature_set", "determined_features"]:
            assert p1.determined_feature_set == p2.determined_feature_set

        elif k == "shap_feature_selector":
            assert p1.shap_feature_selector.__dict__.keys() == p2.shap_feature_selector.__dict__.keys()
            assert p1.shap_feature_selector.shap_top_ns == p2.shap_feature_selector.shap_top_ns
            assert len(p1.shap_feature_selector.shap_values_train_list) == \
                   len(p2.shap_feature_selector.shap_values_train_list)
            for ix, _ in enumerate(p1.shap_feature_selector.shap_values_train_list):
                for a, b in zip(p1.shap_feature_selector.shap_values_train_list[ix],
                                p2.shap_feature_selector.shap_values_train_list[ix]):
                    for aa, bb in zip(a, b):
                        aa: np.ndarray
                        bb: np.ndarray
                        assert (aa == bb).all()

            assert p1.shap_feature_selector.n_features_to_selected_features_list == \
                   p2.shap_feature_selector.n_features_to_selected_features_list

            assert p1.shap_feature_selector.n_features_to_aggregated_features == \
                   p2.shap_feature_selector.n_features_to_aggregated_features

        elif k == "eval_metrics":
            if not (p1.eval_metrics.scoring_metrics_data is None and p2.eval_metrics.scoring_metrics_data is None):
                assert p1.eval_metrics.scoring_metrics_data.equals(
                    p2.eval_metrics.scoring_metrics_data)
                assert p1.eval_metrics.scoring_metrics_data_melted.equals(
                    p2.eval_metrics.scoring_metrics_data_melted)

            assert len(p1.eval_metrics.Xs_benchmark_feature_names_to_dataframes_list) == \
                   len(p2.eval_metrics.Xs_benchmark_feature_names_to_dataframes_list)
            for ix, _ in enumerate(p1.eval_metrics.Xs_benchmark_feature_names_to_dataframes_list):
                assert p1.eval_metrics.Xs_benchmark_feature_names_to_dataframes_list[ix].keys() == \
                       p2.eval_metrics.Xs_benchmark_feature_names_to_dataframes_list[ix].keys()
                for k_item in p1.eval_metrics.Xs_benchmark_feature_names_to_dataframes_list[ix].keys():
                    assert p1.eval_metrics.Xs_benchmark_feature_names_to_dataframes_list[ix][k_item].equals(
                        p2.eval_metrics.Xs_benchmark_feature_names_to_dataframes_list[ix][k_item])

        elif k == "fine_tuner":
            for i in range(len(p1.fine_tuner.classifiers_attributes_data)):
                assert str(p1.fine_tuner.classifiers_attributes_data.iloc[i, 0]) == \
                       str(p2.fine_tuner.classifiers_attributes_data.iloc[i, 0])
                assert str(p1.fine_tuner.classifiers_attributes_data.iloc[i, 1]) == \
                       str(p2.fine_tuner.classifiers_attributes_data.iloc[i, 1])
                assert p1.fine_tuner.classifiers_attributes_data.iloc[i, 2] == \
                       p2.fine_tuner.classifiers_attributes_data.iloc[i, 2], i

            assert all([str(a) == str(b) for a, b in zip(p1.fine_tuner.randomized_search_objects,
                                                         p2.fine_tuner.randomized_search_objects)])

            assert all([str(a) == str(b) for a, b in zip(p1.fine_tuner.best_estimators,
                                                         p2.fine_tuner.best_estimators)])

            assert p1.fine_tuner.param_grid == p2.fine_tuner.param_grid

        elif k in ["ensemble_voting_classifier", "predictions"]:
            # Ensure that ensemble_voting_classifier attribute of both objects are not initialized yet.
            assert (p1.ensemble_voting_classifier is None and p2.ensemble_voting_classifier is None)
            # Ensure that predictions attribute of both objects are not initialized yet.
            assert (p1.predictions is None and p2.predictions is None)

        else:
            raise AttributeError(f"Attribute {k} not handled.")

    log.info(f"All assertions passed! Given Predator objects are the same.")
