import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedShuffleSplit

from ..labels import ClassLabels


def label_proportions(data_param):
    return data_param["Mutation_Effect_Label"].value_counts() / len(data_param)


def get_label_proportions_data(data_prepared, strat_valid_set, random_valid_set):
    compare_props = pd.DataFrame({
        "Overall": label_proportions(data_prepared),
        "Stratified": label_proportions(strat_valid_set),
        "Random": label_proportions(random_valid_set)
    }).sort_index()
    compare_props[f"Rand. %error"] = 100 * compare_props["Random"] / compare_props["Overall"] - 100
    compare_props[f"Strat. %error"] = 100 * compare_props["Stratified"] / compare_props["Overall"] - 100
    compare_props.rename(index={ClassLabels.DISRUPTING: 'Disruptive', ClassLabels.NONDISRUPTING: 'Increasing + No Effect'}, inplace=True)
    return compare_props


# FIXME: prepare *train* data for machine learning
def prepare_data_machine_learning(data, random_seed):
    """

    :param data:
    :param random_seed:
    :return:
    """

    data = data.drop(["UniProt_ID", "Mutation", "Interactor_UniProt_ID"], axis='columns')

    # Shuffle the data
    data_prepared = data.sample(frac=1, random_state=random_seed).reset_index(drop=True).copy()

    # Train and Validation variables
    random_train_set, random_valid_set = train_test_split(
        data_prepared, test_size=0.2, random_state=random_seed)

    split = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=random_seed)
    for train_index, valid_index in split.split(data_prepared, data_prepared["Mutation_Effect_Label"]):
        strat_train_set = data_prepared.iloc[train_index]
        strat_valid_set = data_prepared.iloc[valid_index]

    label_proportions_data = get_label_proportions_data(data_prepared, strat_valid_set, random_valid_set)

    # Declare `X_train`, `y_train`, `X_valid`, `y_valid`
    # All data, i.e. data_prepared
    X = data_prepared.drop(["Mutation_Effect_Label"], axis="columns")
    y = data_prepared["Mutation_Effect_Label"].copy()

    # Stratified version
    X_train = strat_train_set.drop(["Mutation_Effect_Label"], axis="columns")
    y_train = strat_train_set["Mutation_Effect_Label"].copy()
    X_valid = strat_valid_set.drop(["Mutation_Effect_Label"], axis="columns")
    y_valid = strat_valid_set["Mutation_Effect_Label"].copy()

    # Randomized version
    X_train_random = random_train_set.drop(["Mutation_Effect_Label"], axis="columns")
    y_train_random = random_train_set["Mutation_Effect_Label"].copy()
    X_valid_random = random_valid_set.drop(["Mutation_Effect_Label"], axis="columns")
    y_valid_random = random_valid_set["Mutation_Effect_Label"].copy()

    data_materials = {
        "data_prepared": data_prepared,
        "label_proportions_data": label_proportions_data,
        "X": X,
        "y": y,
        "X_train": X_train,
        "y_train": y_train,
        "X_valid": X_valid,
        "y_valid": y_valid,
        "X_train_random": X_train_random,
        "y_train_random": y_train_random,
        "X_valid_random": X_valid_random,
        "y_valid_random": y_valid_random
    }

    return data_materials
