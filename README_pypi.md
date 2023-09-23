[//]: # (# CodonU)

[![PyPI - License](https://img.shields.io/pypi/l/CodonU)](https://opensource.org/licenses/MIT)
[![Downloads](https://static.pepy.tech/personalized-badge/codonu?period=total&units=international_system&left_color=grey&right_color=blue&left_text=Downloads)](https://pepy.tech/project/codonu)
[![CodeFactor](https://www.codefactor.io/repository/github/souradiptoc/codonu/badge/master)](https://www.codefactor.io/repository/github/souradiptoc/codonu/overview/master)
[![Docs](https://img.shields.io/badge/docs-passing-brightgreen)](https://souradiptoc.github.io/CodonU/)

# Introduction

Welcome to CodonU!

This is an integrated package for codon usage analysis.

# Functionalities

Various functionalities can be found below. For gene/genome analysis, this package can:

- For Nucleotide sequences
    - Calculate RSCU
    - Calculate CAI
    - Calculate CBI
    - Calculate ENc
    - Calculate tAI
- For Protein sequences
    - Calculate Aromaticity
    - Calculate gravy

One can also can calculate the multivariate analysis, popularly known as correspondence analysis (COA) for the codons
easily.
Supported calculations are:

- For Nucleotide sequences
    - COA using codon frequency with scale set to length of gene
    - COA using codon RSCU values with scale set to length of gene
- For Protein sequences
    - COA using amino acid frequency with scale set to gravy score
    - COA using amino acid frequency with scale set to aromaticity score

Phylogenetic analysis and tree building now can be done.

Detailed instructions on how to use the functions can be found in
the [examples](https://github.com/SouradiptoC/CodonU/tree/master/Examples)

# Installation

    pip install CodonU

# Future Plans

None. If you would like recommend one, please mail
at [sourochaudhuri@gmail.com](mailto:sourochaudhuri@gmail.com)

