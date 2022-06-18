from datetime import datetime

from ..mylogger import get_handler
import logging


import os
import pathlib

from IPython.display import display

import pandas as pd

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


def save_auc_scores(
        prediction_file_path,
        file_name,
        auc_scores,
        overwrite
):
    folder_path = pathlib.Path(prediction_file_path).parent
    file_name = os.path.join(folder_path, file_name)
    file_date = datetime.today().strftime('%Y-%m-%d')
    file_name = f'{file_name}_{file_date}.csv'

    # Ensure the file is not exists before creating to prevent overwriting.
    if os.path.isfile(file_name) and not overwrite:
        log.error(f"File {file_name} is already exist.\n"
                  "To overwrite existing file, use `overwrite=True`.")

    else:
        auc_scores_data = pd.DataFrame(auc_scores)

        auc_scores_data.index.name = "Method"

        display(auc_scores_data)

        auc_scores_data.to_csv(file_name)

        log.info(f"AUC scores are saved into file {file_name}")


        # data_entries = []
        # with open(file_name, 'w') as file:
        #     for analysis_type, auc_scores in auc_scores.items():
        #         log.debug(f"{analysis_type}: {auc_scores}")
        #         file.write(f"{analysis_type}: {auc_scores}\n")
        #         data_entries.append((analysis_type, auc_scores))
        #
        # log.info(f"AUC scores are saved into file {file_name}")
        #
        # return data_entries

def save_to_excel(prediction_file_path, preliminary_data, file_name, export_flag=True):

    if export_flag:
        folder_path = pathlib.Path(prediction_file_path).parent
        file_name = os.path.join(folder_path, file_name)
        output_file_extension = 'xlsx'
        output_file_date = datetime.today().strftime('%Y-%m-%d')
        output_file_name = f'{file_name}_{output_file_date}.{output_file_extension}'

        # Ensure the file is not exists before creating to prevent overwriting.
        if os.path.isfile(output_file_name):
            log.warning(f'File {output_file_name} is already exist.')

        else:
            # Export
            preliminary_data.to_excel(output_file_name, index=False)
            log.debug(f'{output_file_name} is exported.')

        contents = []
        for column_name in preliminary_data.columns:
            column_name = column_name.replace('/', '_div_')  # replace forward slash with '_div_' if exists.
            filename = "../data/column_descriptions/" + column_name + ".txt"
            with open(filename, 'r') as file:
                content = file.read()  # .split("\n")
                contents.append(content)

        column_description_data = pd.DataFrame({
            "COLUMN_NAME": list(preliminary_data.columns),
            "DESCRIPTION": contents}
        )

        output_file_name = f'{file_name}_{output_file_date}_descriptions.{output_file_extension}'
        save_excel_sheet(output_file_name, column_description_data)
        log.debug(f"descriptions_{output_file_name} is exported.")

    else:
        print('export flag is False')


def save_excel_sheet(excel_path_param, data):
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(excel_path_param, engine='xlsxwriter')

    # Convert the dataframe to an XlsxWriter Excel object.
    data.to_excel(writer, sheet_name='COLUMN_DESCRIPTIONS', index=False)

    workbook = writer.book
    worksheet = writer.sheets['COLUMN_DESCRIPTIONS']

    workbook_format = workbook.add_format({'text_wrap': True})

    # Setting the format
    worksheet.set_column('A:A', 40)
    worksheet.set_column('B:B', 100, workbook_format)

    # Close the Pandas Excel writer and output the Excel file.
    writer.save()
#     workbook.close()
