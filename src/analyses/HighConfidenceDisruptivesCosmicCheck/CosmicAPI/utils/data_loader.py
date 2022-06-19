import os
import pandas as pd
import os.path as op
import glob

from MyAPI import MyCosmicAPI

from utils.log_script import ColorHandler, MyLog
import logging

log = MyLog("debug_runner", level=logging.INFO)
log.addHandler(ColorHandler())

pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)


class DataLoaderHelper:
    def __init__(self, high_confidence_disruptive_data_path):
        self.high_confidence_disruptive_data_path = high_confidence_disruptive_data_path
        self.confidences = self.get_confidences()
        pass

    def get_confidences(self):
        confidences = os.listdir(self.high_confidence_disruptive_data_path)
        return confidences

    def _get_tcga_files(self, confidence):
        tcga_files_folder_path_name = op.join(
            self.high_confidence_disruptive_data_path,
            confidence,
            "*.csv"
        )
        tcga_files = glob.glob(tcga_files_folder_path_name)

        return tcga_files

    @staticmethod
    def _get_scrap_inputs(data_path):
        """
        Returns a tuple: list of input genes and list of input mutations.
        :param data_path:
        :return:
        """
        tcga_data = pd.read_csv(data_path)
        input_genes = tcga_data["GENE"].tolist()
        input_mutations = tcga_data["Mutation"].tolist()

        inputs = {
            "input_genes": input_genes,
            "input_mutations": input_mutations,
        }

        return inputs

    def scrap(self):
        for confidence in self.confidences:
            log.info(f"= = = Confidence: {confidence} = = =")
            tcga_files = self._get_tcga_files(confidence)

            for tcga_file in tcga_files:
                log.info(f"tcga file: {tcga_file}")
                inputs = self._get_scrap_inputs(tcga_file)

                for input_gene, input_mutation in zip(inputs["input_genes"], inputs["input_mutations"]):
                    log.debug(f"RUNNING INPUT GENE: {input_gene}, INPUT MUTATION: {input_mutation}")
                    if str(input_gene) == "nan":
                        log.warning(f"RUNNING INPUT GENE: {input_gene}, INPUT MUTATION: {input_mutation}")
                        log.warning("Skipping NAN value in GENE.")
                        continue

                    MyCosmicAPI(gene=input_gene, mutation=input_mutation)
