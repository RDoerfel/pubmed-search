#%%
from Bio import Entrez, Medline
import pandas as pd
import numpy as np
from datetime import datetime
import os

# %%
def search(query,**kwargs):
    print("Searching papers.")
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(**kwargs,term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
    print("Fetching details.")

    ids = ','.join(id_list)
    Entrez.email = 'your.email@example.com'
    handle = Entrez.efetch(db="pubmed", 
                           id=ids, 
                           rettype="medline",
                           retmode="text")
    results = Medline.parse(handle)
    return results

def store_search(df,result_dir):
    print('Storing results.')
    timestemp = datetime.now().strftime("%m%d%Y_%H%M%S")
    result_name = 'pumbed_{}.xlsx'.format(timestemp)
    resultpath = os.path.join(result_dir,result_name)
    df.to_excel(resultpath,index=False)

def extract_info(papers):
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

def exclude_by_type(dfPapers,exclude_file):
    dfExclusion = pd.read_excel(exclude_file)
    to_exclude = dfExclusion['ExludedPubType'].to_list()
    n_papers = len(dfPapers)
    lKeep = np.zeros([n_papers,1],dtype=bool)
    for i,set_type in enumerate(dfPapers['Type']):
        if(set(set_type) & set(to_exclude) == set()): 
            lKeep[i] = True

    print('Excluding {} papers.'.format(n_papers-lKeep.sum()))
    return dfPapers[lKeep]

def pubmed_search(term,search_params,exclude_file,result_dir):
    results = search(term,**search_params)
    id_list = results['IdList']
    papers = fetch_details(id_list)
    papers = list(papers)
    print('Found {} papers.'.format(len(papers)))

    dfPapers = extract_info(papers)
    dfPapersCleaned = exclude_by_type(dfPapers,exclude_file)
    store_search(dfPapersCleaned,result_dir)

# %%

# %%
