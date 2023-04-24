[//]: # (# CodonU)

![CodonU](https://github.com/SouradiptoC/CodonU/blob/master/images/CODON_U_Background.png)

[![PyPI - License](https://img.shields.io/pypi/l/CodonU)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7857924.svg)](https://doi.org/10.5281/zenodo.7857924)
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
      month        = apr,
      year         = 2023,
      publisher    = {Zenodo},
      version      = {v1.0.2},
      doi          = {10.5281/zenodo.7857924},
      url          = {https://doi.org/10.5281/zenodo.7857924}
    }

# Reference

* J. L. Bennetzen and B. D. Hall, “Codon selection in yeast.,” Journal of Biological Chemistry, vol. 257, no. 6, pp.
  3026–3031, Mar. 1982,
  doi: [10.1016/s0021-9258(19)81068-210.1016/s0021-9258(19)81068-2](https://doi.org/10.1016/S0021-9258(19)81068-2).
* F. Wright, “The ‘effective number of codons’ used in a gene,” Gene, vol. 87, no. 1, pp. 23–29, Mar. 1990,
  doi: [10.1016/0378-1119(90)90491-9](https://doi.org/10.1016/0378-1119(90)90491-9).
* A. Fuglsang, “The ‘effective number of codons’ revisited,” Biochemical and Biophysical Research Communications, vol.
  317, no. 3, pp. 957–964, May 2004, doi: [10.1016/j.bbrc.2004.03.138](https://doi.org/10.1016/j.bbrc.2004.03.138).
* P. M. Sharp and W.-H. Li, “The codon adaptation index-a measure of directional synonymous codon usage bias, and its
  potential applications,” Nucleic Acids Research, vol. 15, no. 3, pp. 1281–1295, 1987,
  doi: [10.1093/nar/15.3.1281](https://doi.org/10.1093/nar/15.3.1281).
* J. Kyte and R. F. Doolittle, “A simple method for displaying the hydropathic character of a protein,” Journal of
  Molecular Biology, vol. 157, no. 1, pp. 105–132, May 1982,
  doi: [10.1016/0022-2836(82)90515-0](https://doi.org/10.1016/0022-2836(82)90515-0).
* J. R. Lobry and C. Gautier, “Hydrophobicity, expressivity and aromaticity are the major trends of amino-acid usage in
  999 Escherichia coli chromosome-encoded genes,” Nucleic Acids Research, vol. 22, no. 15, pp. 3174–3180, Aug. 1994,
  doi: [10.1093/nar/22.15.3174](https://doi.org/10.1093/nar/22.15.3174).
* J. R. Lobry, Multivariate Analyses of Codon Usage Biases. ISTE Press - Elsevier, 2018.
  doi: [10.1016/C2018-0-02165-9](https://doi.org/10.1016/C2018-0-02165-9).
* F. Sievers and D. G. Higgins, “Clustal Omega for making accurate alignments of many protein sequences,” Protein
  Science, vol. 27, no. 1, pp. 135–145, Jan. 2018, doi: [10.1002/pro.3290](https://doi.org/10.1002/pro.3290).
* M. A. Larkin et al., “Clustal W and Clustal X version 2.0,” Bioinformatics, vol. 23, no. 21, pp. 2947–2948, Nov. 2007,
  doi: [10.1093/bioinformatics/btm404](https://doi.org/10.1093/bioinformatics/btm404).
