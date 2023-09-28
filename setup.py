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
VERSION = '1.1.0'
DESC = 'An integrated package for codon usage analysis'

setup(
    name='CodonU',
    version=VERSION,
    description=DESC,
    long_description=read("README_pypi.md"),
    long_description_content_type="text/markdown",
    author=AUTHOR_NAME,
    author_email=AUTHOR_EMAIL,
    url='https://github.com/SouradiptoC/codon_usage',
    project_urls={
        "Bug Tracker": "https://github.com/SouradiptoC/codon_usage/issues",
        "Documentation": "https://codonu.readthedocs.io/en/latest/?badge=latest",
        "Source Code": "https://github.com/SouradiptoC/codon_usage"
    },

    license='MIT License',
    packages=find_packages(exclude=["tests", ".github"]),
    install_requires=["biopython~=1.79", "numpy~=1.24.2", "pandas~=2.1.1",
                      "matplotlib>=3.6.3", "scipy~=1.10.0",
                      "prince~=0.12.1", "gtAI~=1.0.4"],
    keywords=['bioinformatics', 'bioinformatics-analysis', 'bioinformatics-tool', 'codon-usage', 'codon', 'codon-bias',
              'genomic-analysis', 'genome', 'codonW', 'tAI'],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research"
    ]
)
