from oauth2client.service_account import ServiceAccountCredentials
from decouple import config
import numpy as np
import pandas as pd
import gspread
import re
import os
import string
import random


def jsonkey_key_path():
    """get jsonkey file path from environment

    Returns:
        str: path to the .json credential file
    """
    if os.path.exists(".env"):
        jsonkey_path = config("JSONKEY_PATH")
    
    else:
        print("jsonkey path does not exist in environment, setting path as current directory")
        jsonkey_path = os.getcwd()

    return jsonkey_path


def googlesheet_connection(spreadsheetname, jsonpath=os.getcwd(), worksheet=0):
    """function to connect google sheet

    Args:
        spreadsheetname (str): name of spreadsheet contains sample information
        jsonpath (str): path to .json file for googleAPI connection
        worksheet (int): which worksheet to connect

    Returns:
        connected worksheet
    """

    scope = ["https://spreadsheets.google.com/feeds",
             'https://www.googleapis.com/auth/spreadsheets',
             "https://www.googleapis.com/auth/drive.file",
             "https://www.googleapis.com/auth/drive"]

    credential = ServiceAccountCredentials.from_json_keyfile_name(jsonpath, scope)
    client = gspread.authorize(credential)

    # open spreadsheet
    spreadsheet = client.open(spreadsheetname)
    sheet = spreadsheet.get_worksheet(worksheet)
    return sheet


def generate_random_id(length):
    """generate a random file name
    :param length: length of the file name
    """
    chrs = string.ascii_lowercase
    file_name = "".join(random.choice(chrs) for i in range(length))
    return file_name


def get_qc_stat(path, sample_id, level="sample"):
    """get qc statistics from a finished run

    Args:
    path (str): path to pipeline processed result (after re-structure)
    sample_id (str): sample id, e.g. BUKMAP_20190822A
    level (str): experiment/sample - load experiment/sample - level QC statistics

    Returns:
    list object contains qc statistics
    """
    qc_path = ""
    if level == "sample":
        qc_path = os.path.join(path, "Sample_output/QCs")
        qc_stat = pd.read_csv(os.path.join(qc_path, sample_id + ".qc.txt"), sep="\t")
        qc_record = [sample_id, "finshed", "None", path] + qc_stat.iloc[:, 1].tolist()
        return qc_record

    elif level == "experiment":
        qc_path = os.path.join(path, "Experiment_output/QCs")
        # list all qc files in folder
        qc_files = os.listdir(qc_path)
        dfs = []
        for file in qc_files:
            experiment_id = re.sub(".qc.txt", "", file)
            qc_stat = pd.read_csv(os.path.join(qc_path, file), sep="\t")

            # set column names and index
            qc_stat.set_index(list(qc_stat.columns[[0]]), inplace=True)
            qc_stat.columns = [experiment_id]
            dfs.append(qc_stat)

        qc_stat = pd.concat(dfs, axis=1).T
        # insert experiment information
        qc_stat.insert(0, "experiment_id", qc_stat.index)
        qc_stat.insert(1, "sample_id", sample_id, allow_duplicates=True)
        qc_stat.insert(2, "running_stat", "finished", allow_duplicates=True)
        qc_stat.insert(3, "error_log", "None", allow_duplicates=True)
        qc_stat.insert(4, "report_path", path, allow_duplicates=True)
        return qc_stat.values.tolist()
