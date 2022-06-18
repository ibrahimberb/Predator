# Imports
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import LeaveOneOut
# from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import plot_confusion_matrix
from sklearn.metrics import classification_report

from typing import Union

from pandas import DataFrame, Series

AccuracyScore = float
BalancedAccuracyScore = float


def get_default_classifier(random_state=None) -> RandomForestClassifier:
    if random_state is None:
        random_state = 42
    return RandomForestClassifier(random_state=random_state)


def get_cross_validation_option(cv_option):
    """
    A helper function that returns desired cv_option.

    Returns
    -------
        cv_option_param
        n_jobs_param
    """

    # Options for cross-validation cv= parameter.
    if isinstance(cv_option, int):
        cv_option = cv_option
    elif type(cv_option) in [StratifiedKFold, KFold, LeaveOneOut]:
        cv_option = cv_option
    elif cv_option == "skf_5":
        cv_option = StratifiedKFold(shuffle=True, n_splits=5)
    elif cv_option == "skf_10":
        cv_option = StratifiedKFold(shuffle=True, n_splits=10)
    elif cv_option == "kf_5":
        cv_option = KFold(shuffle=True, n_splits=5)
    elif cv_option == "kf_10":
        cv_option = KFold(shuffle=True, n_splits=10)
    elif cv_option == "loocv":
        cv_option = LeaveOneOut()
    else:
        raise ValueError("cv_option value error!")

    return cv_option


def evaluate_cross_val(X_train, y_train, cv_option, n_jobs=-1) -> None:
    """

    :param X_train:
    :param y_train:
    :param cv_option:
    :param n_jobs:
    :return:
    """
    # Cross Validation options
    cv_option = get_cross_validation_option(cv_option)

    # Model
    classifier = get_default_classifier(random_state=42)

    # Cross-validation Balanced Accuracy Scores
    balan_acc_scores = cross_val_score(
        classifier,
        X_train,
        y_train,
        cv=cv_option,
        scoring="balanced_accuracy",
        n_jobs=n_jobs,
    )

    # Cross-validation Accuracy Scores
    acc_scores = cross_val_score(
        classifier, X_train, y_train, cv=cv_option, scoring="accuracy", n_jobs=n_jobs
    )

    # Print scores and averages
    print("Balanced accuracy score AVG : {:.4f}".format(balan_acc_scores.mean()))
    print("Accuracy score AVG          : {:.4f}".format(acc_scores.mean()))


def evaluate_valid(
        clf: RandomForestClassifier,
        X_train: DataFrame,
        y_train: Series,
        X_valid: DataFrame,
        y_valid: Series,
        show_confusion_matrix=False,
) -> Union[AccuracyScore, BalancedAccuracyScore]:
    """
    Evaluates the performance of the given model using  balanced accuracy and accuracy scoring metrics.
    Model is fitted with X_train, y_train and predictions obtained for X_valid. These predictions then
    compared with y_valid.
    :param clf: The model to be used in prediction task.
    :param X_train:
    :param y_train:
    :param X_valid:
    :param y_valid:
    :param show_confusion_matrix:
    :return:
    """
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_valid)

    acc_score = accuracy_score(y_valid, y_pred)
    balanced_acc_score = balanced_accuracy_score(y_valid, y_pred)

    print("Accuracy score\t\t: {:.4f}".format(acc_score))
    print("Balanced accuracy score : {:.4f}".format(balanced_acc_score))

    if show_confusion_matrix:
        plot_confusion_matrix(clf, X_valid, y_valid)

    return acc_score, balanced_acc_score


# def final_scoring(X_train, y_train, X_valid, y_valid) -> None:
#     # Final scoring comparison: X_train_selected_exhaustive_5, y_train with prediction of X_valid_selected_exhaustive_5
#     clf = RandomForestClassifier(random_state=42)
#     clf.fit(X_train, y_train)
#     y_pred = clf.predict(X_valid)
#
#     print("Balanced accuracy score : {:.4f}".format(balanced_accuracy_score(y_valid, y_pred)))
#     print("Accuracy score\t\t: {:.4f}".format(accuracy_score(y_valid, y_pred)))
#
#     plot_confusion_matrix(clf, X_valid, y_valid)


def cross_val_confusion_matrix_via(
        model_param, X_train_param, y_train_param, print_report=False
):
    skf = StratifiedKFold(shuffle=True, n_splits=10)
    y_pred_temp = cross_val_predict(model_param, X_train_param, y_train_param, cv=skf)

    # TODO: dictionary will be a better choice, such as 
    #  d={ClassLabels.DISRUPTING: 'Disrupting', ClassLabels.NONDISRUPTING: 'NoEffect+Increasing'}
    #  so that we don't need to worry about correct label being assigned.
    # label_names = ["Disrupting", "NoEffect+Increasing"]
    label_names = ["LABEL_X", "LABEL_Y"]

    sns.heatmap(
        confusion_matrix(y_train_param, y_pred_temp),
        annot=True,
        fmt="d",
        xticklabels=label_names,
        yticklabels=label_names,
    )
    plt.title(r"$\mathbf{Confusion\ Matrix}$", fontsize=16, fontweight="bold")
    plt.ylabel("Actual", fontsize=16, fontweight="bold")
    plt.xlabel("Predicted", fontsize=16, fontweight="bold")
    plt.show()

    if print_report:
        print(
            classification_report(y_train_param, y_pred_temp, target_names=label_names)
        )
