#%%
from src import io
from src import search

import argparse
import os

#%%
parser = argparse.ArgumentParser(description='Execute Pubmed Search.')
parser.add_argument("-f", required=True, help="Path to json file containing search terms and settings")
 
#%%
if __name__ == "__main__":
    file = 'C:/Git/pubmed-search/data/search_brainage.json'
else:
    args = parser.parse_args() 
    file = args.f

print("Reading: {}".format(file))
search_params = io.read_json(file)
term = search_params['Terms']
print("Terms: {}".format(term))

settings = search_params['Settings']
print("Settings: {}".format(settings))

exclude_file = search_params['Exclude']

result_dir = search_params['Result']

search.pubmed_search(term,settings,exclude_file,result_dir)

# %%
