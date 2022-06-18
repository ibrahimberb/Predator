import pandas as pd
from pandas import DataFrame
from IPython.display import display
import seaborn as sns

from ..mylogger import get_handler
import logging

from ..labels import ClassLabels

handler = get_handler()
handler_simple = get_handler(log_type="simple")

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)

log_simple = logging.getLogger(__name__)
log_simple.handlers[:] = []
log_simple.addHandler(handler_simple)
log_simple.setLevel(logging.DEBUG)


def display_data(data: DataFrame, head_row=3):
    log_simple.debug(f"Data dimensions: {data.shape}")
    display(data.head(head_row))


def display_label_counts(data_param):
    """
    Display a dataframe that contains label categories and their counts.
    """
    label_counts = pd.DataFrame(data_param["Mutation_Effect_Label"].value_counts())
    label_counts.reset_index(inplace=True)
    label_counts.columns = ["Mutation_Effect_Label", "Counts"]
    if len(label_counts[label_counts == ClassLabels.DISRUPTING]) > len(label_counts[label_counts == ClassLabels.NONDISRUPTING]):
        # this is index of the data, so the order will be like this.
        if ClassLabels.DISRUPTING == 0:
            raise NotImplementedError

        else:
            raise NotImplementedError

    else: # there are more nondisruptive entry
        if ClassLabels.DISRUPTING == 0:
            raise NotImplementedError

        else:
            label_counts.rename(
                {
                    0: 'Disrupting',
                    1: 'Increasing + No Effect',
                },
                inplace=True
            )

    display(label_counts)


# TODO: unittest
def display_labels(data_param):
    """
    Display a dataframe that contains label categories.
    """
    label_counts = pd.DataFrame(data_param["Mutation_Effect_Label"].value_counts().index)
    label_counts.columns = ["Mutation_Effect_Label"]
    display(label_counts)


def visualize_label_counts(data_param, label_name_param="Mutation_Effect_Label"):
    # noinspection SpellCheckingInspection
    sns.set(style="white", font_scale=1.15)  # white, dark, whitegrid, darkgrid, ticks
    val_counts = data_param[label_name_param].value_counts().sort_index()
    val_counts = val_counts.rename({ClassLabels.DISRUPTING: 'Disrupting', ClassLabels.NONDISRUPTING: 'Increasing + No Effect'})
    log.debug(f"Label counts:\n{val_counts}")
    ax = sns.barplot(x=val_counts.index,
                     y=val_counts,
                     palette="ch:s=-.2,r=.6")
    ax.set_title('Disrupting vs Increasing & No Effect')  # ch:s=-.2,r=.6, ocean
    ax.set_ylabel('Value counts')
