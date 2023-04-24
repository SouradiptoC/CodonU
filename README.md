[//]: # (# CodonU)

![CodonU](https://github.com/SouradiptoC/CodonU/blob/master/images/CODON_U_Background.png)

[![PyPI - License](https://img.shields.io/pypi/l/CodonU)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/536583655.svg)](https://zenodo.org/badge/latestdoi/536583655)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/CodonU)](https://pypi.org/project/CodonU)
[![PyPI](https://img.shields.io/pypi/v/CodonU)](https://pypi.org/project/CodonU)
[![PyPI - Format](https://img.shields.io/pypi/format/CodonU)](https://pypi.org/project/CodonU)
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

# Graphics!!

Also, can generate beautiful graphics for publication purpose or otherwise. Some plots are:

ENc plot for human chromosome 2

![ENc plot for human chromosome 2](https://github.com/SouradiptoC/CodonU/blob/master/images/ENc_plot_Human%20Cr%202.png)

Neutrality plot for human chromosome 2

![Neutrality plot for human chromosome 2](https://github.com/SouradiptoC/CodonU/blob/master/images/Neutrality_plot_Human%20Cr%202.png)

Correspondence analysis of protein frequency using gravy score

![Correspondence analysis of protein frequency using gravy score](https://github.com/SouradiptoC/CodonU/blob/master/images/Multivariate_analysis_aa_gravy_agnetis.png)

Examples for other plots cans be found in the [images](https://github.com/SouradiptoC/CodonU/tree/master/images)

# Installation

    pip install CodonU

# Future Plans

None. If you would like recommend one, please mail
at [sourochaudhuri@gmail.com](mailto:sourochaudhuri@gmail.com)

# Citation

Please cite the work as

    @software{souradipto_choudhuri_2023_7657797,
      author       = {Souradipto Choudhuri},
      title        = {CodonU},
      month        = feb,
      year         = 2023,
      publisher    = {Zenodo},
      version      = {v1.0.0},
      doi          = {10.5281/zenodo.7657797},
      url          = {https://doi.org/10.5281/zenodo.7657797}
    }

