#%%
import json
import shutil
from yattag import Doc
import pandas as pd

# %%

def df_to_html(data):
    """Convert data to html
    
    Parameters
    ----------
    data : pandas data frame
        Data to convert to html.
    
    Returns
    -------
    str
        Html code.
    """
    doc, tag, text, line = Doc().ttl()
    with tag('head'):
        with tag('meta', name="viewport", content="width=device-width, initial-scale=1"):
            pass
        with tag('link', rel="stylesheet", href="html/style.css"):
            pass

    with tag('script', type="text/javascript", src="html/collapsible.js"):
        pass

    with tag('body'):
        for i,row in data.iterrows():
            title = row['Titles']
            abstract = row['Abstracts']
            authors = row['Authors']
            pmid = row['PMID']
            date = row['Date']
            journal = row['Journal']
            with tag('button', type="button", klass='collapsible'):
                with tag('h3', id = 'heading'):
                    text(f'{title}')
            with tag('div', klass='content'):
                with tag('p'):
                    with tag('b'):
                        text(f'Journal: ')
                    text(f'{journal}')
                with tag('p'):
                    with tag('b'):
                        text(f'Authors: ')
                    text(f'{authors}')
                with tag('p'):
                    with tag('b'):
                        text(f'Date: ')
                    text(f'{date}')
                with tag('p'):
                    with tag('b'):
                        text(f'PMID: ')
                    text(f'{pmid}')
                with tag('p'):
                    text(f'{authors}')
                with tag('p'):
                    text(f'{abstract}')
        with tag('script'):
            text('makeCollabsible()')    
    result = doc.getvalue()
    return result

# %%
