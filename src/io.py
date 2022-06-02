import json
import pandas as pd

def read_json(file):
    with open(file) as json_file:
        data = json.load(json_file)
    return data

def read_excel_column(file,key):
    dfExclusion = pd.read_excel(file)
    column = dfExclusion[key].to_list()
    return column
