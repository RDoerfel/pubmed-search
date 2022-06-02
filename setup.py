from setuptools import find_packages, setup

setup(
    name="pmsearch",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
                'pmsearch=pmsearch.pubmed_search:run',
        ]
    },
    install_requires=[
        'biopython==1.78',
        'pandas==1.3.5',
        'numpy==1.21.2',
        'openpyxl==3.0.9'
    ]
)
