import json
import shutil
import pandas as pd
from datetime import datetime
import os

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

def store_search(df,result_dir):
    """Store dataframe with time stemp

    Parameters
    ----------
    df : DataFrame
        Results to store. 
    result_dir : str
        diretory to store results
    """
    print('Storing results.')
    timestemp = datetime.now().strftime("%m%d%Y_%H%M%S")
    result_name = 'pubmed_{}.xlsx'.format(timestemp)
    resultpath = os.path.join(result_dir,result_name)
    df.to_excel(resultpath,index=False)

def store_html(html,result_dir):
    """Store dataframe with time stemp

    Parameters
    ----------
    html : string
        html string to write to file. 
    result_dir : str
        diretory to store results
    """
    print('Storing results.')
    with open(os.path.join(result_dir,'index.html'), 'wb') as f:
        f.write(bytes(html, encoding='utf8'))
        html_path = os.path.join('html/')
        dest = os.path.join(result_dir,'html')
        shutil.copytree(html_path, dest)
