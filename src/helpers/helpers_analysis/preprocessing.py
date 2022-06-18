from pandas import DataFrame
from ..mylogger import get_handler
import logging

handler = get_handler('module')

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


def process_snv(data) -> DataFrame:
    """
    Given a SNV dataframe, do the following processes and return it.

        1) Filter for `Variant_Classification`:
            Keep `Missense_Mutation` only.

        2) Filter for `SWISSPROT` column:
            Drop nan (empty) values.

        3) Filter for `SWISSPROT` column:
            Drop nan (empty) values.

        4) Processing `HGVSp_Short` column:
            For mutation, exclude 'p' in the beginning of string.

    Parameters
    ----------
        data : <DataFrame>
            SNV dataframe to be processed.

    Returns
    -------
        data : <DataFrame>
            Processed version of input dataframe.

    """

    # 1) Filter for `Variant_Classification` column, and keep `Missense_Mutation` only.
    data = data[data["Variant_Classification"] == "Missense_Mutation"]

    # 2) Filter for 'SWISSPROT' column and drop nan (empty) values.
    data = data.dropna(subset=["SWISSPROT"])

    # 3) Filter for 'HGVSp_Short' column and drop nan (empty) values.
    data = data.dropna(subset=["HGVSp_Short"])

    # 4) Processing `HGVSp_Short` column, removing preceeding 'p.'.
    data["HGVSp_Short"] = data["HGVSp_Short"].apply(lambda x: str(x).replace('p.', ''))

    # Reset index of the dataframe to avoid any possible errors
    data.reset_index(drop=True, inplace=True)

    return data


def shorten_patient_ids(data) -> DataFrame:
    """
    Shorten the Patient IDs by taking first 12 char.
    Example shortened id: TCGA-D8-A1XY
    """

    data['Tumor_Sample_Barcode'] = data['Tumor_Sample_Barcode'].apply(lambda x: x[:12])

    return data


def simplify_snv_data(data) -> DataFrame:
    """
    Given SNV data, simplify it into four columns: "Hugo_Symbol", "SWISSPROT", "HGVSp_Short" and "Tumor_Sample_Barcode".
    """

    data = data[["Hugo_Symbol", "SWISSPROT", "HGVSp_Short", "Tumor_Sample_Barcode"]].copy(deep=True)

    return data


def simplify_elaspic_data(data) -> DataFrame:
    """
    Given ELASPIC results data, simplify it into three columns: "UniProt_ID", "Mutation", "Interactor_UniProt_ID".
    """

    data = data[["UniProt_ID", "Mutation", "Interactor_UniProt_ID"]].copy(deep=True)

    return data


def remove_duplicated_entries(data) -> DataFrame:
    """
    Removed duplicated rows in given data.
    """

    log.debug(f'{data.duplicated().sum()} duplicated entries are removed.')
    data = data.drop_duplicates(keep="first")

    return data


def convert_primary_isomer(column_name: str, data: DataFrame) -> DataFrame:
    """
    Converts proteins into primary form representation (dash-free from) in given column name of the given dataframe.
    E.g.
        P16473-2 → P16473

    Parameters
    ----------
        column_name : <string>
            Name of the column where protein is stored.

        data : <DataFrame>
            The dataframe whose proteins will be processed in `column_name` column.

    Returns
    -------
        data : <DataFrame>
            Processed version of input dataframe.

    """

    # Protein names will be converted dashed-free version, if they contain.
    data[column_name] = data[column_name].apply(lambda x: x.split('-')[0])

    return data
