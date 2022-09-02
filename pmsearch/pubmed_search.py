#%%
from pmsearch.to_html import df_to_html
from . import version
from . import io
from . import search
import argparse
import shutil
# %%
def pubmed_search(terms,settings,result_dir,exclude_file=None,include=None,html=False):
    """Run pubmed search for specified settings

    Parameters
    ----------
    terms : str
        Terms for pubmed search. Could be "(Pubmed) AND (Search)".
    settings : dict
        Settings for pubmed search.
    result_dir : str
        Directory to store results in
    exclude_file : str, optional
        Filename with journal types to exclude, by default None
    include : str, optional
        String with comma separeted PMIDs to look for in results, by default None
    html : bool, optional
        Apply html formatting to results, by default False
    """
    print("=================================Settings=================================")

    print("Terms: {}".format(terms))    
    print("Settings: {}".format(settings))
    
    if(exclude_file):
        exclude = io.read_excel_column(exclude_file,'ExludedPubType')
        print("Excluding: {}".format(exclude) )
    else: 
        exclude = []
        print("Excluding: {}".format("None") )

    if(include):
        include = include.split(',')
        print("Check for: {}".format(include))
    else:
        include = []
        print("Check for: {}".format("None"))

    print("=================================Run=================================")
    id_list = search.search(terms,**settings)
    papers = search.fetch_details(id_list)
    print('Found {} papers.'.format(len(papers)))

    dfPapers = search.extract_info(papers)
    dfPapersCleaned = search.exclude_by_type(dfPapers,exclude)
    search.check_for_papers(dfPapersCleaned,include)
    io.store_search(dfPapersCleaned,result_dir)

    if(html):
        print('Storing results as html.')
        html_text = df_to_html(dfPapersCleaned)
        io.store_html(html_text,result_dir)

def run():
    """Run command."""
    version()
    parser = argparse.ArgumentParser(description='Execute PubMed search.')

    parser.add_argument("-f", required=True, help="Path to json file containing search terms and settings")
    parser.add_argument("--html", required=False, default=False, help="Apply html formatting to results")

    args = parser.parse_args() 
    file = args.f

    search_params = io.read_json(file)
    pubmed_search(**search_params)

if __name__ == '__main__':
    run()
