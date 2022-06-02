import json
import pandas as pd

def read_json(file):
    """Read Json file

    Parameters
    ----------
    file : str
        json file.

    Returns
    -------
    dict
        Content of json file
    """
    with open(file) as json_file:
        data = json.load(json_file)
    return data

def read_excel_column(file,key):
    """Read a column from an excel file

    Parameters
    ----------
    file : str
        The Excel file
    key : str
        the column to read

    Returns
    -------
    list
        The column.
    """
    dfExclusion = pd.read_excel(file)
    column = dfExclusion[key].to_list()
    return column
