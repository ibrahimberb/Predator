import os
import pandas as pd
from matplotlib import pyplot as plt
from datetime import datetime
import os.path as op

import sys
sys.path.insert(0, "../../../../../src/")

from helpers.labels import ClassLabels

# os.system("python Mutation_pred.py -i 1a22.pdb -c A -r 182 -w C -m A -d 0")
PREDATOR_BENCHMARK_DIR = "."
TRAIN_DATA_PATH = op.join(PREDATOR_BENCHMARK_DIR, "../../../common/train_data_with_PDB.csv")


def save_prediction_data(predator_benchmark_dir, prediction_file_name, prediction_data):
    file_date = datetime.now().strftime("%Y-%m-%d")
    prediction_file_name = "{}_{}.csv".format(prediction_file_name, file_date)
    prediction_data.to_csv(op.join(predator_benchmark_dir, prediction_file_name), index=False)
    print("Prediction data `{}`is exported.".format(op.join(predator_benchmark_dir, prediction_file_name)))


def save_errors(errors, file_path):
    with open(file_path, "w") as f:
        for error in errors:
            f.write("Error type: %s\n" % error)
            
    print("Errors written to file: " + file_path)


def read_train_data(file_name):
    df = pd.read_csv(file_name)
    return df


def run_SAAMBE_3D(pdb_file, chain, residue_position, wildtype_aa, mutant_aa, model_type):
    system_command = (
            "python "
            + "Mutation_pred.py "
            + "-i {} ".format(pdb_file)
            + "-c {} ".format(chain)
            + "-r {} ".format(residue_position)
            + "-w {} ".format(wildtype_aa)
            + "-m {} ".format(mutant_aa)
            + "-d {} ".format(model_type)
    )

    print("Running system command: {}".format(system_command))
    output = os.popen(system_command).read().strip()

    return output


def main():
    # result = run_SAAMBE_3D(pdb_file="1a22.pdb", chain="A", residue_position="182", wildtype_aa="C", mutant_aa="A", model_type="0")
    # print("Result: {}".format(result))

    output_label_mapping = {
        "Disruptive": ClassLabels.DISRUPTING,
        "Nondisruptive": ClassLabels.NONDISRUPTING,
    }

    train_data = read_train_data(TRAIN_DATA_PATH)
    prediction_data = train_data.copy()

    errors = []
    error_messages = set()
    
    counts = {
        "errors": 0,
        "true_predictions": 0,
        "false_predictions": 0,
    }

    predictions = []

    for index, row in train_data.iterrows():
        print(" Processing row {} of {} ".format(index, train_data.shape[0]).center(100, "-"))

        pdb_file = row["Template_cath_id_pdb"]
        pdb_path = op.join(PREDATOR_BENCHMARK_DIR, "pdb_files", pdb_file + ".pdb")

        protein = row["UniProt_ID"]
        interactor = row["Interactor_UniProt_ID"]
        mutation = row["Mutation"]
        chain = row["Chain_id"]
        wildtype_aa = mutation[0]
        mutant_aa = mutation[-1]
        residue_position = mutation[1:-1]
        model_type = 0
        actual_class_label = row["Mutation_Effect_Label"]

        print("pdb_file: {}".format(pdb_file))
        print("protein: {}".format(protein))
        print("mutation: {}".format(mutation))
        print("interactor: {}".format(interactor))

        output = run_SAAMBE_3D(pdb_path, chain, residue_position, wildtype_aa, mutant_aa, model_type)

        if output not in ["Nondisruptive", "Disruptive"]:
            error_string = ("\033[91m> Error: {} {} {} {} {} {}\n {} \033[00m".format(
                pdb_file, chain, residue_position, wildtype_aa, mutant_aa, model_type, output)
            )
            print(error_string)
            errors.append(error_string)
            error_messages.add(output)
            counts["errors"] += 1
            predictions.append(None)  # We don't have a prediction for this row
            continue

        output_class_label = output_label_mapping[output]

        if output_class_label == actual_class_label:
            print("\033[92mPredicted: {} ({}) Actual: {} \033[00m".format(output, output_class_label, actual_class_label))
            counts["true_predictions"] += 1

        else:
            print("\033[91mPredicted: {} ({}) Actual: {} \033[00m".format(output, output_class_label, actual_class_label))
            counts["false_predictions"] += 1

        # Add prediction to predictions list
        predictions.append(output_class_label)

        # print(f"Predicted: {output} ({output_class_label}) Actual: {actual_class_label}")

    # Save the results
    prediction_data["Prediction"] = predictions
    save_prediction_data(PREDATOR_BENCHMARK_DIR, "SAAMBE_3D_predictions", prediction_data)
    save_errors(error_messages, "{}\\errors_messages.txt".format(PREDATOR_BENCHMARK_DIR))

    # plot the distribution of the counts in dictionary "counts"
    plt.title("Prediction Results")
    plt.bar(range(len(counts)), list(counts.values()), align="center")
    plt.xticks(range(len(counts)), list(counts.keys()))
    plt.savefig("{}\\prediction_results.png".format(PREDATOR_BENCHMARK_DIR))
    plt.show()


if __name__ == '__main__':
    main()
