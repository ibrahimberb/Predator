from datetime import datetime


def executed_on():
    print("\033[32m{}\033[0m".format(datetime.now().strftime("%B %d, %Y %H:%M:%S")))
