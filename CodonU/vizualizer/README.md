# Vizualizer

The package is designed for vizualizing the data. This is documentation for vizualizer package.
Functions for general use:

- [plot_enc](#plot_enc)
- [plot_pr2](#plot_pr2)
- [plot_neutrality](#plot_neutrality)

Functions for advanced use:

- [_enc](#_enc)
- [_plot_enc](#_plot_enc)
- [_plot_pr2](#_plot_pr2)
- [_plot_neutrality](#_plot_neutrality)

# Detailed Description of Functions

Detailed description of functions are given below.

## Detailed description for functions of general usage

### `plot_enc`

Plots ENc curve from given fasta file

- `handle`: Handle to the file, or the filename as a string
- `genetic_table_num`: Genetic table number for codon table
- `min_len_threshold`: Minimum length of nucleotide sequence to be considered as gene
- `organism_name`: Name of organism (optional)
- `save_image`: Options for saving the image (optional)
- `folder_path`: Folder path where image should be saved (optional)
- `gene_analysis`: Option if gene analysis (True) or genome analysis (False) (optional)

In the plot, the super title is name of the organism or only name if `organism_name` is not given. The title of the
curve is number of genes or genome. On the right side, a color bar is provided for better understanding.
See [example](#enc-plot).

### `plot_pr2`

Plots A3/AT3 values against G3/GC3 values from given fasta file

- `handle`: Handle to the file, or the filename as a string
- `min_len_threshold`: Minimum length of nucleotide sequence to be considered as gene
- `organism_name`: Name of organism (optional)
- `save_image`: Options for saving the image (optional)
- `folder_path`: Folder path where image should be saved (optional)
- `gene_analysis`: Option if gene analysis (True) or genome analysis (False) (optional)

In the plot, the super title is name of the organism or only name if `organism_name` is not given. The title of the
curve is number of genes or genome. Each quadrant has the count for the occurrence of purines and pyrimidines. Again the
counts can be seen in the lower left side. See [example](#pr2-plot).

### `plot_neutrality`

Plots neutrality plot from given fasta file

- `handle`: Handle to the file, or the filename as a string
- `min_len_threshold`: Minimum length of nucleotide sequence to be considered as gene
- `organism_name`: Name of organism (optional)
- `save_image`: Options for saving the image (optional)
- `folder_path`: Folder path where image should be saved (optional)
- `gene_analysis`: Option if gene analysis (True) or genome analysis (False) (optional)

In the plot, the super title is name of the organism or only name if `organism_name` is not given. The title of the
curve is number of genes or genome and the R^2 value. On the right side, a color bar is provided for better
understanding. See [example](#neutrality-plot).

## Detailed descriptions of functions for advanced usage

### `_enc`

Calculates Theoretical ENC value based on Wright (1989) ([click here](https://doi.org/10.1016/0378-1119(90)90491-9)).

- `x`: GC3 value

Returns the ENc value.

### `_plot_enc`

Plots ENc value against GC3 values

- `enc_val_lst`: Values of ENc (in range 0 to 1)
- `gc_val_lst`: Values of GC3 (in range 0 to 1)
- `organism_name`: Name of organism (optional)
- `save_image`: Options for saving the image (optional)
- `folder_path`: Folder path where image should be saved (optional)
- `gene_analysis`: Option if gene analysis (True) or genome analysis (False) (optional)

### `_plot_pr2`

Plots A3/AT3 values against G3/GC3 values

- `gc_val_lst`: Values of GC3 (in range 0 to 1)
- `at_val_lst`: Values of AT3 (in range 0 to 1)
- `g3_val_lst`: Values of G3 (in range 0 to 1)
- `a3_val_lst`: Values of A3 (in range 0 to 1)
- `organism_name`: Name of organism (optional)
- `save_image`: Options for saving the image (optional)
- `folder_path`: Folder path where image should be saved (optional)
- `gene_analysis`: Option if gene analysis (True) or genome analysis (False) (optional)

### `_plot_neutrality`

Plots the neutrality plot

- `gc12_lst`: The list containing values of G or C at 1 or 2 positions
- `gc3_lst`: The list containing values of G or C ar 3 positions
- `organism_name`: Name of organism (optional)
- `save_image`: Options for saving the image (optional)
- `folder_path`: Folder path where image should be saved (optional)
- `gene_analysis`: Option if gene analysis (True) or genome analysis (False) (optional)

## Examples

### ENc Plot

#### ENc plot for _Staphylococcus agnetis_

![ENc plot for Staphylococcus agnetis](https://github.com/SouradiptoC/CodonU/blob/master/CodonU/images/ENc_plot_Staphylococcus%20agnetis.png)

#### ENc plot for Human chromosome 2

![ENc plot for human chromosome 2](https://github.com/SouradiptoC/CodonU/blob/master/CodonU/images/ENc_plot_Human%20Cr%202.png)

### PR2 Plot

#### PR2 plot for _Staphylococcus agnetis_

![PR2 plot for Staphylococcus agnetis](https://github.com/SouradiptoC/CodonU/blob/master/CodonU/images/PR2_plot_Staphylococcus%20agnetis.png)

#### PR2 plot for Human chromosome 2

![PR2 plot for human chromosome 2](https://github.com/SouradiptoC/CodonU/blob/master/CodonU/images/PR2_plot_Human%20Cr%202.png)

### Neutrality Plot

The color-bar depends on the span of GC3 or GC12. As can be seen below, the color-bar for _agnetis_ is based on the span
of GC3 whereas for _argenteus_ it is based on span of GC12.

#### Neutrality plot for _Staphylococcus agnetis_

![Neutrality plot for Staphylococcus agnetis](https://github.com/SouradiptoC/CodonU/blob/master/CodonU/images/Neutrality_plot_Staphylococcus%20agnetis.png)

#### Neutrality plot for _Staphylococcus argenteus_

![Neutrality plot for Staphylococcus argenteus](https://github.com/SouradiptoC/CodonU/blob/master/CodonU/images/Neutrality_plot_Staphylococcus%20argenteus.png)

#### Neutrality plot for Human chromosome 2

![Neutrality plot for human chromosome 2](https://github.com/SouradiptoC/CodonU/blob/master/CodonU/images/Neutrality_plot_Human%20Cr%202.png)

<p align="center">&copy; 2023 Souradipto Choudhuri</p>