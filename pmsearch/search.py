#%%
from Bio import Entrez, Medline
import pandas as pd
import numpy as np

# %%
def search(query,**kwargs):
    """Search pubmed with queries

    Parameters
    ----------
    query : str
        Search terms.

    Returns
    -------
    list
        PMIDs for matching papers
    """
    print("Searching papers.")
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(**kwargs,term=query)
    results = Entrez.read(handle)
    return results['IdList']

def fetch_details(id_list):
    """Fetch details for provided PMIDs

    Parameters
    ----------
    id_list : list
        PMIDs to fetch details for

    Returns
    -------
    list
        List with details for identified papers
    """
    print("Fetching details.")

    ids = ','.join(id_list)
    Entrez.email = 'your.email@example.com'
    handle = Entrez.efetch(db="pubmed", 
                           id=ids, 
                           rettype="medline",
                           retmode="text")
    results = Medline.parse(handle)
    return list(results)

def extract_info(papers):
    """Extract information from list of papers and store in dataframe. 

    Parameters
    ----------
    papers : list
        Papers to extract information from.

    Returns
    -------
    DataFrame
        Extracted information.
    """
    lTitles = []
    lAbstracts = []
    lAuthors = []
    lDate = []
    lPMID = []
    lJournal = []
    lType = []

    for i, paper in enumerate(papers):
        lAuthors.append(paper.get('FAU'))
        lTitles.append(paper.get('TI'))
        lAbstracts.append(paper.get('AB'))
        lDate.append(paper.get('DP'))
        lPMID.append(paper.get('PMID'))
        lJournal.append(paper.get('JT'))
        lType.append(paper.get('PT'))
        entries = list(zip(lAuthors, lTitles, lAbstracts, lDate, lPMID, lJournal, lType))
    names = ['Authors', 'Titles', 'Abstracts', 'Date', 'PMID','Journal','Type']
    dfPapers = pd.DataFrame(entries,columns=names)
    return dfPapers

def exclude_by_type(dfPapers,exclude):
    """Exclude papers from search based on Journal Type

    Parameters
    ----------
    dfPapers : DataFrame
        Papers.
    exclude : list
        Journal types to exclude

    Returns
    -------
    DataFrame
        Cleaned papers.
    """
    n_papers = len(dfPapers)
    lKeep = np.zeros([n_papers,1],dtype=bool)
    for i,set_type in enumerate(dfPapers['Type']):
        if(set(set_type) & set(exclude) == set()): 
            lKeep[i] = True

    print('Excluding {} papers.'.format(n_papers-lKeep.sum()))
    return dfPapers[lKeep]

def check_for_papers(df,include):
    """Check if specific papers are in search result.

    Parameters
    ----------
    df : DataFrame
        Papers.
    include : list
        PMIDs from papers to check for
    """
    same = set(df['PMID']) & set(include)
    print("{} out of {} where found in search.".format(len(same),len(include)))
    print(same)

