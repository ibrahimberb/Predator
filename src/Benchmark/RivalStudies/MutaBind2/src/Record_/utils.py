from datetime import datetime
import os.path as op


def save_prediction_data(benchmark_dir, prediction_file_name, prediction_data):
    file_date = datetime.now().strftime("%Y-%m-%d")
    prediction_file_name = "{}_{}.csv".format(prediction_file_name, file_date)
    prediction_data.to_csv(op.join(benchmark_dir, prediction_file_name), index=False)
    print("Prediction data `{}`is exported.".format(op.join(benchmark_dir, prediction_file_name)))
