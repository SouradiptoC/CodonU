# Extractor

The package is designed for extracting the data. This is documentation for extractor package.
Functions for general use:

- [extract_cds](#extract_cds)
- [extract_cds_lst](#extract_cds_lst)
- [extract_exome](#extract_exome)
- [extract_prot](#extract_prot)

Functions for advanced use:

- [extract_cds_seq](#extract_cds_seq)
- [extract_prot_seq](#extract_prot_seq)

# Detailed Description of Functions

Detailed description of functions are given below.

## Detailed description for functions of general usage

### `extract_cds`

Extracts the CDS as a Sequence Record object from a given Sequence Record object

- `record`: Original Sequence Record object from where the CDS is to be extracted
- `feature_location`: The location of CDS, which is a Feature Location object
- `cds_no`: Number of CDS

For details about SeqRecord object [click here](https://biopython.org/docs/1.75/api/Bio.SeqRecord.html). For details
about FeatureLocation object [click here](https://biopython.org/docs/1.75/api/Bio.SeqFeature.html). The function returns
a new Sequence Record object containing the CDS

### `extract_cds_lst`

Extracts the list of features if their type is CDS

- `record`: Original Sequence Record object from where the CDS is to be extracted

For details about SeqRecord object [click here](https://biopython.org/docs/1.75/api/Bio.SeqRecord.html). Returns a tuple
containing FeatureLocation objects. For details
about FeatureLocation object [click here](https://biopython.org/docs/1.75/api/Bio.SeqFeature.html).

### `extract_exome`

Extracts the exome from given nucleotides

- `nuc_file_path`: The path to the nucleotide file
- `organism_name`: Name of the organism

Returns the exome as a Sequence Record object. For details about SeqRecord
object [click here](https://biopython.org/docs/1.75/api/Bio.SeqRecord.html).

### `extract_prot`

Extracts protein sequences from given sequence feature object.

- `feature`: The CDS, which is a SeqFeature object
- `organism_name`: Name of the organism
- `cds_no`: Number of the CDS

For details about SeqFeature object [click here](https://biopython.org/docs/1.75/api/Bio.SeqFeature.html). The function
returns the protein sequence as a SeqRecord object. For details about SeqRecord
object [click here](https://biopython.org/docs/1.75/api/Bio.SeqRecord.html).

## Detailed descriptions of functions for advanced usage

### `extract_cds_seq`

Extracts the CDS sequence from a given sequence

- `seq`: Provided sequence
- `feature_location`: The location of the CDS

`feature_location` is an object of type FeatureLocation. For details
about FeatureLocation object [click here](https://biopython.org/docs/1.75/api/Bio.SeqFeature.html). The function returns
the sequence of CDS which is of type Seq. For details about Seq
object [click here](https://biopython.org/docs/1.75/api/Bio.Seq.html).

### `extract_prot_seq`

Extracts the protein sequence reported in the report for the provided cds.

- `feature`: The CDS, which is a SeqFeature object

For details about SeqFeature object [click here](https://biopython.org/docs/1.75/api/Bio.SeqFeature.html). The function
returns the protein sequence, which is of type Seq. For details about Seq
object [click here](https://biopython.org/docs/1.75/api/Bio.Seq.html).

<p align="center">&copy; 2023 Souradipto Choudhuri</p>