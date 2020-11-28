from oauth2client.service_account import ServiceAccountCredentials
from pprint import pprint
import numpy as np
import pandas as pd
import gspread
import re
import os


def jsonkey_key_path():
    """get jsonkey file path from environment

    Returns:
        str: path to the .json credential file
    """

    if "JSONKEY_PATH" in os.environ:
        jsonkey_path = os.environ["JSONKEY_PATH"]
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
