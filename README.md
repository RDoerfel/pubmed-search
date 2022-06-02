# pubmed-search

Command line application to search Pubmed. This tool utilizes the API provided by [https://biopython.org/docs/latest/api/Bio.Entrez.html]() to communicate with PubMed.

### Getting started

You can install this package either from source or directly from GitHub. This package was developed using `python==3.9`.

#### Installation

Via GitHub:

```
pip install git+https://github.com/RDoerfel/pubmed-search
```

From Source:

```
git clone https://github.com/RDoerfel/pubmed-search.git
cd pubmed-search
pip install .
```

I recommend using the package in its own onda environment to make sure that everything works.

### Usage

#### Define Search Parameters

The search parameters are defined using a json file, like it is provided in `data/search.json`.

##### Mandatory Settings:

`"terms":"(PubMed) AND (Search)"` - The PubMed search query to execute.

`"settings":{"db":"pubmed"}` - The database to use (pubmed is recommended)

`"result_dir":"<path_to_result_dir>"` - Path to store the results.

##### Optional Settings:

```
"settings":
    {
    "db":"pubmed", 
    "sort":"relevance", 
    "mindate":"2010/01/01",
    "maxdate":"2022/05/20",
    "retmax":1000,
    "retmode":"xml"}
```

Settings as specified in [https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch]().

`"exclude_file":"<path_to_exclude_file>"` - The path to the file containing Journal Types to exclude from search

`"include":"PMID1,PMID2"` - list of PMIDs to check if they are included in results.

### Execute search

Run pubmed search using: 

`pmsearch -f <path_to_search.json>`
