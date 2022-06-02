from pmsearch import __app_name__, __version__
from pmsearch import pubmed_search

def version():
    """Return Version"""
    print("Running {} version {}".format(__app_name__,__version__))

if __name__ == "__main__":
    version()
    pubmed_search.run()
