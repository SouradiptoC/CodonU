# Analyzer

The package is designed for analyzing the data. This is documentation for analyzer package.
Functions for general use:

- [calculate_rscu](#calculate_rscu)
- [calculate_cai](#calculate_cai)
- [calculate_cbi](#calculate_cbi)
- [calculate_enc](#calculate_enc)

Functions for advanced use:

- [g3](#g3)
- [a3](#a3)
- [gc_123](#gc_123)
- [at_123](#at_123)
- [filter_reference](#filter_reference)
- [syn_codons](#syn_codons)
- [sf_vals](#sf_vals)
- [cbi](#cbi)
- [enc](#enc)

# Detailed Description of Functions

Detailed description of functions are given below.

## Detailed description for functions of general usage

### `calculate_rscu`

Calculates rscu values for each codon. It is inspired by the works of Sharp and Li (1987). For
details [click here](https://doi.org/10.1093/nar/15.3.1281).

- `records`: The generator object containing sequence object
- `genetic_code_num`: Genetic table number for codon
  table ([click here](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for details)
- `min_len_threshold`: Minimum length of nucleotide sequence to be considered as gene (Default 200)
- `gene_analysis`: Option if gene analysis (True) or genome analysis (False) (Default false)

The function returns dictionary containing codon and rscu value pairs, if `gene_analysis` is false. Else will return the
dictionary containing gene name and the dictionary which contains codon and rscu value pairs.

### `calculate_cai`

Calculates codon adaptation index values for each codon. It is inspired by the works of Sharp and Li (1987). For
details [click here](https://doi.org/10.1093/nar/15.3.1281).

- `records`: The generator object containing sequence object
- `genetic_code_num`: Genetic table number for codon
  table ([click here](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for details)
- `min_len_threshold`: Minimum length of nucleotide sequence to be considered as gene (Default 200)
- `gene_analysis`: Option if gene analysis (True) or genome analysis (False) (Default false)

The function returns dictionary containing codon and cai value pairs, if `gene_analysis` is false. Else will return the
dictionary containing gene name and the dictionary which contains codon and cai value pairs.

### `calculate_cbi`

Calculates codon bias index values for each amino acid. It is inspired by the works of Bennetzen and Hall (1982). For
details [click here](https://doi.org/10.1016/S0021-9258(19)81068-2).

- `records`: The generator object containing sequence object
- `genetic_code_num`: Genetic table number for codon
  table ([click here](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for details)
- `min_len_threshold`: Minimum length of nucleotide sequence to be considered as gene (Default 200)
- `gene_analysis`: Option if gene analysis (True) or genome analysis (False) (Default false)

The function returns dictionary containing amino acid and cbi value, optimal codon pairs, if `gene_analysis` is false.
Else it will return the dictionary containing gene name and the dictionary which contains amino acid and cbi value,
optimal codon pairs.

### `calculate_enc`

Calculates effective number of codons value for a given sequences. It is inspired by the works of Wright (1989). For
details [click here](https://doi.org/10.1016/0378-1119(90)90491-9). Some revisions to the original method have been done
for certain classes of
protein based on the works of Fuglsang (2004). For details [click here](https://doi.org/10.1016/j.bbrc.2004.03.138).

- `records`: The generator object containing sequence object
- `genetic_code_num`: Genetic table number for codon
  table ([click here](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for details)
- `min_len_threshold`: Minimum length of nucleotide sequence to be considered as gene (Default 200)
- `gene_analysis`: Option if gene analysis (True) or genome analysis (False) (Default false)

The function returns the effective number of codons value for given sequence if `gene_analysis` is false. Else it will
return a dictionary containing gene number and corresponding ENc value

## Detailed descriptions of functions for advanced usage

### `g3`

Calculates percentage of G content for third position

- `seq`: Provided sequence

Returns the percentage of G content for third position

### `a3`

Calculates percentage of A content for third position

- `seq`: Provided sequence

Returns the percentage of A content for third position

### `gc_123`

Calculate G+C content: total, for first, second and third positions

- `seq`: Provided sequence

Returns the G+C percentage for the entire sequence, and the three codon positions

### `at_123`

Calculate A+T content: total, for first, second and third positions

- `seq`: Provided sequence

Returns the A+T percentage for the entire sequence, and the three codon positions

### `filter_reference`

Filters the list of reference based on given threshold of length

- `records`: A generator object holding the sequence objects
- `min_len_threshold`: Minimum length of nucleotide sequence to be considered as gene

Returns the list of usable sequences

### `syn_codons`

Creates the protein, codon dictionary where protein is key

- `codon_table`: The codon table

Returns the dictionary of protein, codon pair having protein as key

### `sf_vals`

Creates the sf value and protein dictionary where sf value is key

- `codon_table`: The codon table

Returns the dictionary of sf value and protein dictionary where sf value is key

### `cbi`

Calculates codon bias index (CBI) for a given protein seq based on Bennetzen and Hall (1982). For
details [click here](https://doi.org/10.1016/S0021-9258(19)81068-2).

- `prot_seq`: The Amino Acid
- `reference`: List of reference nucleotide sequences
- `genetic_code`: Genetic table number for codon
  table ([click here](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for details)

This function may raise the following warnings

- `NoSynonymousCodonWarning`: When there is no synonymous codons
- `MissingCodonWarning`: When no codons translate to provided Amino acid

Returns a tuple of CBI vals and the optimal codons

### `enc`

Calculates Effective number of codons (Enc) for a given list of sequences based on Wright (1989)
([click here](https://doi.org/10.1016/0378-1119(90)90491-9)) and Fuglsang
(2004) ([click here](https://doi.org/10.1016/j.bbrc.2004.03.138)).

- `references`: List of reference nucleotide sequences
- `genetic_code`: Genetic table number for codon table

This function may raise the following warnings and errors

- `MissingCodonWarning`: When no codons translate to provided Amino acid
- `NoProteinError`: If there is no codon for a certain set of amino acid

Returns calculated Enc value for the sequence(s)

## References

- Bennetzen, J.L. and Hall, B.D. (1982) “Codon selection in yeast.,” Journal of Biological Chemistry, 257(6), pp.
  3026–3031. Available at: https://doi.org/10.1016/s0021-9258(19)81068-2
- Fuglsang, A. (2004) “The ‘effective number of
  codons’ revisited,” Biochemical and Biophysical Research Communications, 317(3), pp. 957–964. Available
  at: https://doi.org/10.1016/j.bbrc.2004.03.138
- Sharp, P.M. and Li, W.-H. (1987) “The codon adaptation index-a measure
  of directional synonymous codon usage bias, and its potential applications,” Nucleic Acids Research, 15(3), pp.
  1281–1295. Available at: https://doi.org/10.1093/nar/15.3.1281
- Wright, F. (1990) “The ‘effective number of codons’ used
  in a gene,” Gene, 87(1), pp. 23–29. Available at: https://doi.org/10.1016/0378-1119(90)90491-9.

<p align="center">&copy; 2023 Souradipto Choudhuri</p>