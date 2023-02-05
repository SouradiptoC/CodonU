# Codon Usage Analysis Warnings

This package holds all the warnings used in other srcs. The list of warnings are given below.

- [CodonUsageWarning](#codonusagewarning)
- [ApiWarning](#apiwarning)
- [EmailWarning](#emailwarning)
- [MissingCodonWarning](#missingcodonwarning)
- [NoSynonymousCodonWarning](#nosynonymouscodonwarning)

# Detailed Description

### `CodonUsageWarning`

Accounts for all warnings while executing the files

### `ApiWarning`

Occurs when no API key is provided

### `EmailWarning`

Occurs when no email id is provided

### `MissingCodonWarning`

Occurs when no codon in the given reference sequence list translates to a certain amino acid

- `aa`: The amino acid which is not translated

### `NoSynonymousCodonWarning`

Occurs when only one codon in the given reference sequence list translates to a certain amino acid

- `aa`: The amino acid which is translated from one codon