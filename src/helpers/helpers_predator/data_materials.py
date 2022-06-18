from .loaders import load_train_data, load_tcga_data
from .preprocessing import preprocess_train_data, preprocess_tcga_data

from .prepare_data_for_ML import prepare_data_machine_learning

from tqdm.notebook import tqdm

from ..mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class DataMaterials(dict):

    def __init__(self, n_experiment, random_seeds):
        super().__init__()
        self.n_experiment = n_experiment
        self.random_seeds = random_seeds

    def initialize_model_data_materials(self):
        """
        Instantiate datasets into labels and features for machine learning.
        Includes train and validations splits as well.
        """
        log.debug("Initializing model datasets ..")
        self["prepared_dataframes"] = []
        self["label_proportions_dataframes"] = []
        self["Xs"] = []
        self["ys"] = []
        self["Xs_train"] = []
        self["ys_train"] = []
        self["Xs_valid"] = []
        self["ys_valid"] = []
        self["Xs_train_random"] = []
        self["ys_train_random"] = []
        self["Xs_valid_random"] = []
        self["ys_valid_random"] = []

    def append_data_materials(self, data_materials: dict):
        log.debug("Appending data materials ..")
        self["prepared_dataframes"].append(data_materials['data_prepared'])
        self["label_proportions_dataframes"].append(data_materials['label_proportions_data'])
        self["Xs"].append(data_materials['X'])
        self["ys"].append(data_materials['y'])
        self["Xs_train"].append(data_materials['X_train'])
        self["ys_train"].append(data_materials['y_train'])
        self["Xs_valid"].append(data_materials['X_valid'])
        self["ys_valid"].append(data_materials['y_valid'])
        self["Xs_train_random"].append(data_materials['X_train_random'])
        self["ys_train_random"].append(data_materials['y_train_random'])
        self["Xs_valid_random"].append(data_materials['X_valid_random'])
        self["ys_valid_random"].append(data_materials['y_valid_random'])

    def initialize_feature_selected_data_materials(self, n_top: int):
        log.debug(f"Initialize feature selected data materials [n_top={n_top}]..")
        self[f"Xs_shap_{n_top}"] = []
        self[f"Xs_train_shap_{n_top}"] = []
        self[f"Xs_valid_shap_{n_top}"] = []
        # Single feature columns
        self[f"Xs_provean"] = []
        self[f"Xs_train_provean"] = []
        self[f"Xs_ddG"] = []
        self[f"Xs_train_ddG"] = []

    def initialize_train_dataset(self, project_common_file_dir, initial_columns_path, mutations_path):
        log.debug("Initialize `train_data` ..")
        log.debug("Initialize `train_data_processed` ..")
        self["train_data"] = load_train_data(project_common_file_dir, mutations_path)
        self["train_data_processed"] = preprocess_train_data(
            self["train_data"], initial_columns_path
        )

    def initialize_target_datasets(self, project_common_file_dir, initial_columns_path, tcga_code_path_pairs):
        for tcga_code_path_pair in tcga_code_path_pairs:
            tcga_name, tcga_path = tcga_code_path_pair
            tcga_name = tcga_name.lower()
            log.debug(f"Initialize `{tcga_name}` ..")
            log.debug(f"Initialize `target_{tcga_name}_data` ..")
            self[f"{tcga_name}"] = load_tcga_data(project_common_file_dir, tcga_path)
            self[f"target_{tcga_name}_data"] = preprocess_tcga_data(
                self[f"{tcga_name}"], initial_columns_path
            )

    def prepare_model_data_materials(self):
        log.debug("preparing datasets for ML ..")
        self.initialize_model_data_materials()
        for i in tqdm(range(self.n_experiment)):
            data_materials = prepare_data_machine_learning(
                self["sampled_train_data_list"][i], random_seed=self.random_seeds[i]
            )
            self.append_data_materials(data_materials)

    def append_feature_selected_data_materials(self, n_top, selected_features):
        for exp in range(self.n_experiment):
            # SHAP
            self[f"Xs_shap_{n_top}"].append(self[f"Xs"][exp][selected_features])
            self[f"Xs_train_shap_{n_top}"].append(self[f"Xs_train"][exp][selected_features])
            self[f"Xs_valid_shap_{n_top}"].append(self[f"Xs_valid"][exp][selected_features])
            # single feature columns
            self[f"Xs_provean"].append(self["Xs"][exp][["Provean_score"]])
            self[f"Xs_train_provean"].append(self["Xs_train"][exp][["Provean_score"]])
            self[f"Xs_ddG"].append(self["Xs"][exp][["Final_ddG"]])
            self[f"Xs_train_ddG"].append(self["Xs_train"][exp][["Final_ddG"]])

    def initialize_target_data_materials(self, determined_features, tcga_code_path_pairs):
        log.debug("Initializing target data materials ..")
        log.debug(f"Determined features: {determined_features}")
        for tcga_code_path_pair in tcga_code_path_pairs:
            tcga, tcga_path = tcga_code_path_pair
            tcga = tcga.lower()
            log.debug(f"Declaring Xs_{tcga} data materials  ..")
            Xs_tcga = []
            for _ in range(self.n_experiment):
                Xs_tcga.append(
                    self[f"target_{tcga}_data"][determined_features].copy(deep=True)
                )
            self[f"Xs_{tcga}"] = Xs_tcga
