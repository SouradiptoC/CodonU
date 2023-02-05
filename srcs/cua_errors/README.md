# Codon Usage Analysis Errors

This package holds all the warnings used in other srcs. The list of errors are given below.

- [CodonUsageError](#codonusageerror)
- [FileNotEmptyError](#filenotemptyerror)
- [internal_stop_codon_err](#internalstopcodonerror)
- [no_protein_err](#noproteinerror)
- [nucleotide_err](#nucleotideerror)
- [ter_seq_err](#terseqerror)

# Detailed Description

### `CodonUsageError`

Accounts for all error that can happen during executing the files

### `FileNotEmptyError`

Occurs when a given file to write is not empty

- `file_name`: The path or name of the file

### `InternalStopCodonError`

Occurs when internal stop codon is present in the genome

- `position`: Position of codon (begging from 0)

### `NoProteinError`

Occurs when a complete category of amino acid based on sf values is not translated by the provided sequence

- `seq`: The sequence which is incompatible

### `NucleotideError`

Occurs when an ambiguous or invalid nucleotide is present in genome

- `code`: error code (1 for ambiguous, 2 for invalid)
- `position`: Position of codon (begging from 0)

### `TerSeqError`

Occurs when terminal codons are not valid

- `code`: error code (1 for starting, -1 for ending)