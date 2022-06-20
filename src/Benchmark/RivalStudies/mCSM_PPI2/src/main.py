from MyAPI import MyMCSMPP2API
from src.utils.errors import MissionAbortError
from src.utils.resources import terminate_firefox_processes
from utils.log_script import ColorHandler
from utils.record import TrainDataIndexes
import pandas as pd

import logging

log = logging.Logger("debug_runner", level=logging.DEBUG)
log.addHandler(ColorHandler())

# input_mutation_effect_label = 1
# input_protein = "P01116"
# input_mutation = "Y40C"
# input_interactor = "P50749"
# input_pdb_id = "3DDC"
# input_chain_id = "A"


def main():
    TRAIN_DATA_WITH_PDB_PATH = "../../../common/train_data_with_PDB.csv"
    train_data = pd.read_csv(TRAIN_DATA_WITH_PDB_PATH)
    train_data_indexes = TrainDataIndexes()

    for index, row in train_data.iterrows():

        log.debug(f"  Running index {index}  ".center(100, "-"))

        if index in train_data_indexes.successful_indexes:
            continue

        if index in train_data_indexes.failed_indexes:
            continue

        input_mutation_effect_label = row["Mutation_Effect_Label"]  # 1
        input_protein = row["UniProt_ID"]  # "P01116"
        input_mutation = row["Mutation"]  # "Y40C"
        input_interactor = row["Interactor_UniProt_ID"]  # "P50749"
        input_pdb_id = row["Template_cath_id_pdb"]  # "3DDC"
        input_chain_id = row["Chain_id"]  # "A"

        try:
            MyMCSMPP2API(
                mutation_effect_label=input_mutation_effect_label,
                protein=input_protein,
                mutation=input_mutation,
                interactor=input_interactor,
                pdb_id=input_pdb_id,
                chain_id=input_chain_id,
                upload_our_pdb=True,
            )

            train_data_indexes.add_successful_index(index)

        except MissionAbortError:
            train_data_indexes.add_failed_index(index)
            continue


if __name__ == '__main__':
    main()
    # terminate_firefox_processes()

# entry_1 = Entry(
#     mutation_effect_label=1,
#     protein="protein_1",
#     mutation="mutation_1",
#     interactor="interactor_1",
#     pdb_id="pdb_id_1",
#     chain_id="chain_id_1",
#     submitted=1
# )
#
# entry_2 = Entry(
#     mutation_effect_label=2,
#     protein="protein_2",
#     mutation="mutation_2",
#     interactor="interactor_2",
#     pdb_id="pdb_id_2",
#     chain_id="chain_id_2",
#     submitted=2
# )


# record = Record()
#
#
# print(record.record_data.shape)
# print(record.record_data)
#
# print("adding entry 1")
# record.add_entry(entry_1)
# print(record.record_data.shape)
# print(record.record_data)
#
# print("adding entry 2")
# record.add_entry(entry_2)
# print(record.record_data.shape)
# print(record.record_data)
#
