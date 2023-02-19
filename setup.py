import io
import os
from setuptools import setup, find_packages


def read(*paths, **kwargs):
    """
    Read the contents of a text file safely.

    """
    content = ""
    with io.open(
            os.path.join(os.path.dirname(__file__), *paths),
            encoding=kwargs.get("encoding", "utf8"),
    ) as open_file:
        content = open_file.read().strip()
    return content


def read_requirements(path):
    return [
        line.strip()
        for line in read(path).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


PACK_NAME = 'CodonU'
AUTHOR_NAME = 'Souradipto Choudhuri'
AUTHOR_EMAIL = 'sourochaudhuri@gmail.com'
VERSION = '1.0.0'
DESC = 'This package is designed for helping in genomic analysis'
LONG_DESC = 'The package can easily can help in various ways from fetching genebank files from NCBI with the help of ' \
            'accession id to calculating cai, cbi, rscu, enc values and also can generate good quality graphics ' \
            'such as enc plot, or correspondence analysis plot'

setup(
    name='CodonU',
    version=VERSION,
    url='https://github.com/SouradiptoC/codon_usage',
    author=AUTHOR_NAME,
    author_email=AUTHOR_EMAIL,
    description=DESC,
    long_description=LONG_DESC,
    packages=find_packages(exclude=["tests", ".github"]),
    install_requirements=read_requirements('requirements.txt'),
    keywords=['bioinformatics', 'bioinformatics-analysis', 'bioinformatics-tool', 'codon-usage', 'codon', 'codon-bias',
              'genomic-analysis', 'genome', 'souradipto choudhuri'],
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Bioinformaticians",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows"
    ]
)
