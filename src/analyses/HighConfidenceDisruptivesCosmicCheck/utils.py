from datetime import datetime
from pathlib import Path
import os.path as op
from typing import Union

import numpy as np
from pandas import DataFrame
import pandas as pd

from IPython.display import display

from src.analyses.HighConfidenceDisruptivesCosmicCheck.CosmicAPI.utils.misc import get_residue_position
from src.helpers.helpers_analysis.gene_id_retrieval import GeneIDFetcher

from src.helpers.helpers_analysis.convert_primary_isomer import convert_primary_isomer

from src.helpers.labels import ClassLabels

from src.helpers.mylogger import get_handler
import logging

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.WARNING)

log.propagate = False

UNIPROT_GENE_MAPPING_PATH = "../../helpers/helpers_analysis/gene_retrieval/UNIPROT_GENE_MAPPING.csv"

# set pandas max column names unlimited
pd.set_option('display.max_columns', None)

# set width of the output
pd.set_option('display.width', None)


class HighConfidenceDisruptiveMutationsHelper:
    def __init__(self, tcga: str, data: DataFrame, confidence: float):
        """
        :param tcga: The TCGA cohort name
        :param data: Prediction data
        :param confidence: Disruptive confidence should be between 0.50 and 1.
        """
        if confidence < 0.5:
            raise ValueError("If confidence is less then 0.50, it is not predicted as disruptive.")

        self.tcga = tcga.upper()
        self.data = data
        self.confidence = confidence
        self.gene_id_fetcher = GeneIDFetcher(UNIPROT_GENE_MAPPING_PATH)

        self.high_confidence_data = None
        self.prepare_high_confidence_disruptive_mutations()

    def get_high_confidence_disruptive_mutations(self):
        return self.high_confidence_data

    def prepare_high_confidence_disruptive_mutations(self):
        """
        Confidence is between 0 and 1
        """

        # Since Median probability is the probability of being class 1, we need to look what that `1` stands for.
        if ClassLabels.DISRUPTING == 0:
            self.data["Disruptive_probability"] = 1 - self.data["Median_Probability"]

        # The median probability means the disruptive class!
        elif ClassLabels.DISRUPTING == 1:
            self.data["Disruptive_probability"] = self.data["Median_Probability"]

        else:
            raise ValueError

        high_confidence_data = self.data[self.data["Disruptive_probability"] >= self.confidence]

        assert len(high_confidence_data["Prediction"].unique()) == 1, high_confidence_data["Prediction"].nunique()
        assert high_confidence_data["Prediction"].unique()[0] == ClassLabels.DISRUPTING, high_confidence_data["Prediction"].unique()[0]

        high_confidence_data.insert(
            loc=0,
            column="TCGA",
            value=self.tcga
        )

        high_confidence_data.insert(
            loc=1,
            column="GENE",
            value=high_confidence_data["UniProt_ID"].apply(lambda x: self.gene_id_fetcher.fetch(x))
        )

        high_confidence_data.insert(
            loc=4,
            column="INTERACTOR_GENE",
            value=high_confidence_data["Interactor_UniProt_ID"].apply(lambda x: self.gene_id_fetcher.fetch(x))
        )

        self.high_confidence_data = high_confidence_data

    def extract_high_confidence_disruptive_mutations(self, view=False):
        folder_path = op.join(
            "HighConfidenceDisruptiveData", f"confidence_{self.confidence:.2f}"
        )
        file_date = datetime.today().strftime('%Y-%m-%d')
        Path(f"{folder_path}").mkdir(parents=True, exist_ok=True)
        file_name = f"{self.tcga}_confidence_{self.confidence:.2f}_{file_date}.csv"
        file_path = op.join(folder_path, file_name)

        if op.isfile(file_path):
            raise FileExistsError(f"You already have the file {file_path}")

        high_confidence_data = self.get_high_confidence_disruptive_mutations()

        if view:
            display(high_confidence_data)

        high_confidence_data.to_csv(file_path, index=False)

        print(f"{self.tcga} data is extracted to {file_path} successfully.")


class CosmicResultsAttaching:
    def __init__(self, cosmic_results_data):
        self.cosmic_results_data = cosmic_results_data

    @staticmethod
    def find_in_cosmic_results(
            cosmic_results_data: DataFrame,
            gene: str,
            mut: str
    ) -> dict:

        query = cosmic_results_data[
            (cosmic_results_data["GENE"] == gene) &
            (cosmic_results_data["RESIDUE_POSITION"] == int(get_residue_position(mut)))
            ]

        if query.empty:
            query_result = {
                "CGC_status": "NOT_FOUND",
                "most_significant_codon_tier": "NOT_FOUND",
            }

        else:
            [CGC_status] = query["CGC_STATUS"]
            [most_significant_codon_tier] = query["MOST_SIGNIFICANT_CODON_TIER"]

            query_result = {
                "CGC_status": CGC_status,
                "most_significant_codon_tier": most_significant_codon_tier,
            }

        return query_result

    def attach_results(
            self,
            tcga_prediction_data: DataFrame,
    ) -> DataFrame:

        tcga_prediction_data_cosmic_results = tcga_prediction_data.copy()
        cosmic_results = tcga_prediction_data_cosmic_results.apply(
            lambda row: self.find_in_cosmic_results(
                cosmic_results_data=self.cosmic_results_data,
                gene=row["GENE"],
                mut=row["Mutation"]
            ), axis=1
        )

        tcga_prediction_data_cosmic_results["CGC_status"] = cosmic_results.apply(
            lambda row: row["CGC_status"]
        )

        tcga_prediction_data_cosmic_results["MOST_SIGNIFICANT_CODON_TIER"] = cosmic_results.apply(
            lambda row: row["most_significant_codon_tier"]
        )

        return tcga_prediction_data_cosmic_results


class ProveanScoreAttaching:
    """
    It retrieves from downloaded Interface datasets (Merged_Results).
    """

    def __init__(
            self,
            tcga_cosmic_results_data: DataFrame,
            tcga: str,
            tcga_elaspic_results_data_path: Union[str, Path]
    ):
        self.tcga = tcga.lower()
        assert self.tcga == tcga_cosmic_results_data["TCGA"].unique()[0].lower()
        self.tcga_cosmic_results_data = tcga_cosmic_results_data
        self.tcga_data_with_features = pd.read_csv(tcga_elaspic_results_data_path, sep="\t")
        self.tcga_data_with_features = convert_primary_isomer(
            column_name="Interactor_UniProt_ID", data=self.tcga_data_with_features
        )
        self.provean_attached_data = None

    def attach_provean_scores(self, provean_loc=9):
        retrieved_provean_scores = []
        for index, row in self.tcga_cosmic_results_data.iterrows():
            retrieved_provean_scores.append(
                self._get_provean_score(
                    protein=row["UniProt_ID"],
                    mutation=row["Mutation"],
                    interactor=row["Interactor_UniProt_ID"]
                )
            )

        provean_attached_data = self.tcga_cosmic_results_data.copy()

        provean_attached_data.insert(
            loc=provean_loc,
            column="PROVEAN",
            value=retrieved_provean_scores,
        )

        self.provean_attached_data = provean_attached_data

        return self.provean_attached_data

    def _get_provean_score(self, protein, mutation, interactor):
        """
        I know Provean score is not related with interactor,
        but here I will find correct entry.
        """
        query = self.tcga_data_with_features[
            (self.tcga_data_with_features["UniProt_ID"] == protein) &
            (self.tcga_data_with_features["Mutation"] == mutation) &
            (self.tcga_data_with_features["Interactor_UniProt_ID"] == interactor)
            ].copy()

        try:
            query["Provean_score"] = query["Provean_score"].astype(float)
        except ValueError:
            contains_only_none = list(
                query["Provean_score"].values
            ).count("None") == len(query["Provean_score"])

            if contains_only_none:
                log.error(
                    f"\n{protein=}"
                    f"\n{mutation=}"
                    f"\n{interactor=}"
                    f"\nIt contains None only."
                )
                # display(query)
                return "N/A"

            else:
                raise

        if len(query) != 1:
            triplets_columns = ["UniProt_ID", "Mutation", "Interactor_UniProt_ID"]
            if query[triplets_columns].duplicated(keep=False).all():
                # no problem
                try:
                    [provean_score] = query["Provean_score"].unique()
                except ValueError:
                    # Take average of provean scores.
                    log.warning(
                        f"\n{protein=}"
                        f"\n{mutation=}"
                        f"\n{interactor=}"
                        f"\nTook average of {query['Provean_score'].unique()}"
                    )
                    provean_score = np.mean(query["Provean_score"].unique())

            else:
                # print(f"{protein=}")
                # print(f"{mutation=}")
                # print(f"{interactor=}")
                # display(query[["UniProt_ID", "Mutation", "Interactor_UniProt_ID", "Provean_score"]])
                raise WillHandleLaterError

        else:
            [provean_score] = query["Provean_score"]

        return provean_score


class WillHandleLaterError(Exception):
    """I will worry about it later."""
