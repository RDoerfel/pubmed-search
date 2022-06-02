#%%
from src import io
from src import search

import argparse

# %%
def pubmed_search(terms,settings,result_dir,exclude_file=None,include=None):
    print("Terms: {}".format(terms))    
    print("Settings: {}".format(settings))
    
    if(exclude_file):
        exclude = io.read_excel_column(exclude_file,'ExludedPubType')
        print("Excluding: {}".format(exclude) )
    else: 
        exclude = []
        print("Excluding: {}".format("None") )

    if(include):
        print("Check for: {}".format(include))
        include = include.split(',')
    else:
        include = []
        print("Check for: {}".format("None"))

    results = search.search(terms,**settings)
    id_list = results['IdList']
    papers = search.fetch_details(id_list)
    papers = list(papers)
    print('Found {} papers.'.format(len(papers)))

    dfPapers = search.extract_info(papers)
    dfPapersCleaned = search.exclude_by_type(dfPapers,exclude)
    search.check_for_papers(dfPapersCleaned,include)
    search.store_search(dfPapersCleaned,result_dir)

def run():
    parser = argparse.ArgumentParser(description='Execute Pubmed Search.')

    parser.add_argument("-f", required=True, help="Path to json file containing search terms and settings")

    args = parser.parse_args() 
    file = args.f
    print(file)
    search_params = io.read_json(file)
    pubmed_search(**search_params)
    
# %%
if __name__ == "__main__":
    run()
