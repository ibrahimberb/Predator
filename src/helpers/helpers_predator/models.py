from .machine_learning_utils import get_default_classifier
from ..mylogger import get_handler
import logging
from tqdm.notebook import tqdm
from sklearn.base import clone

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class DefaultModels(list):
    def __init__(self, n_experiment):
        super().__init__()
        for _ in range(n_experiment):
            self.append(get_default_classifier(random_state=42))


class TunedModels(list):
    def __init__(self, best_estimators):
        super().__init__()
        for estimator in best_estimators:
            self.append(estimator)


class QualifiedModels(list):
    def __init__(self, qualified_estimators):
        super().__init__()
        assert qualified_estimators is not None
        for estimator in qualified_estimators:
            self.append(estimator)


class FinalizedModels(list):
    def __init__(self, tuned_models):
        super().__init__()
        for estimator in tuned_models:
            estimator_cloned = clone(estimator)
            self.append(estimator_cloned)

    def fit_all(self, data_materials, determined_feature_set):
        for exp, estimator in tqdm(enumerate(self), total=len(self)):
            estimator.fit(data_materials[f"Xs_{determined_feature_set}"][exp],
                          data_materials[f"ys"][exp])


class EnsembleVotingClassifier:
    def __init__(self, models, voting):
        log.debug("Initializing EnsambledVotingClassifier.")
        log.debug(f"Voting mode: {voting}")
        self.models = models
        self.voting = voting

    def predict_hard_voting(self, Xs):
        tcga_predictions = []
        for estimator, X in zip(self.models, Xs):
            log.debug(f"Current estimator: {estimator}")
            log.debug(f"X shape: {X.shape}")
            prediction = estimator.predict(X)
            tcga_predictions.append(prediction)

        return tcga_predictions

    def predict_soft_voting(self, Xs):
        tcga_predictions_probabilities = []
        for estimator, X in zip(self.models, Xs):
            log.debug(f"Current estimator: {estimator}")
            log.debug(f"X shape: {X.shape}")
            prediction_probability = estimator.predict_proba(X)
            # Ensure classes are in following order so that the first probabilities
            # belonging to class 0 (whatever that is), and second probabilities belonging to class 1 (whatever that it).
            assert list(estimator.classes_) == [0, 1]
            tcga_predictions_probabilities.append(prediction_probability)

        return tcga_predictions_probabilities

    def predict(self, Xs):
        if self.voting == "hard":
            return self.predict_hard_voting(Xs)
        elif self.voting == "soft":
            return self.predict_soft_voting(Xs)
        else:
            raise ValueError(f"Invalid attribute. Expected self.voting values are "
                             f"`hard` and `soft`, got {self.voting}")

