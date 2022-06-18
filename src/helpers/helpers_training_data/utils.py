# from ..mylogger import get_handler
# import logging
#
# import pandas as pd
#
# handler = get_handler(log_type="module")
#
# log = logging.getLogger(__name__)
# log.handlers[:] = []
# log.addHandler(handler)
# log.setLevel(logging.DEBUG)
#
#
# def read_intact_mutation_data(intact_file_path):
#     # Read Original IntAct Mutation Data
#     mutations_data = pd.read_table(intact_file_path, delimiter="\t")
#     log.debug(f"Size of dataframe: {mutations_data.shape}")
#
#
# def filter_homo_sapiens(mutations_data):
#     # Filtering the data where "Affected protein organism" is "9606 - Homo sapiens".
#     mutations_homo_sapiens_data = mutations_data[
#         mutations_data['Affected protein organism'] == "9606 - Homo sapiens"
#     ]
#     log.debug(f"Size of dataframe: {mutations_homo_sapiens_data.shape}")
#
