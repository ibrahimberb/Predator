import urllib.request
import os
from datetime import datetime
import os.path as op


# download the pdb file from Protein Data Bank
def download_pdb(pdb_id):
    # print the status of the download
    print("Downloading " + pdb_id + "...")
    url = "http://www.rcsb.org/pdb/files/" + pdb_id + ".pdb"
    file_name = pdb_id + ".pdb"
    
    # create a folder named pdb_files if it doesn't exist
    if not os.path.exists("pdb_files"):
        os.makedirs("pdb_files")

    file_path = os.path.join("pdb_files", file_name)
        
    if not os.path.exists(file_path):
        urllib.request.urlretrieve(url, file_path)
        
    else:
        print("File already exists")
        return


def save_prediction_data(predator_benchmark_dir, prediction_file_name, prediction_data):
    file_date = datetime.now().strftime("%Y-%m-%d")
    prediction_file_name = "{}_{}.csv".format(prediction_file_name, file_date)
    prediction_data.to_csv(op.join(predator_benchmark_dir, prediction_file_name), index=False)
    print("Prediction data `{}`is exported.".format(op.join(predator_benchmark_dir, prediction_file_name)))


def save_errors(errors, file_path):
    with open(file_path, "w") as f:
        for error in errors:
            f.write("%s\n" % error)
            
    print("Errors written to file: " + file_path)
