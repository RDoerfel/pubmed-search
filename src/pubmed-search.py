#%%
from src import io
from src import search

import argparse
import os

#%%
if __name__ == "__main__":
    search_params = io.read_json('C:/Git/pubmed-search/data/search_brainage.json')
    term = search_params['Terms']
    settings = search_params['Settings']
    exclude_file = search_params['Exclude']
    result_dir = search_params['Result']
    
    search.pubmed_search(term,settings,exclude_file,result_dir)

# %%
