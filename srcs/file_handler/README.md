# Extractor

The package is designed for reading and writing the data. This is documentation for file handler package.
Functions for general use:

- [get_gb](#get_gb)
- [make_dir](#make_dir)
- [read_file](#read_file)
- [set_entrez_param](#set_entrez_param)
- [write_exome_fasta](#write_exome_fasta)
- [write_nucleotide_fasta](#write_nucleotide_fasta)
- [write_protein_fasta](#write_protein_fasta)

# Detailed Description of Functions

Detailed description of functions are given below.

## Detailed description for functions of general usage

### `get_gb`

Retrieves the Sequence Record object from a given accession number

- `accession_id`: Provided accession number

Returns the retrieved genebank (gb) file as a SeqRecord object. For details about SeqRecord
object [click here](https://biopython.org/docs/1.75/api/Bio.SeqRecord.html).

### `make_dir`

Makes a directory if not present already

- `path`: Path of the directory

### `read_file`

Returns a dataframe from given csv file

- `file_name`: Name or path to csv file

Returns the read cvs file in the pandas DataFrame object. For details about
pd.DataFrame [click here](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html).

### `set_entrez_param`

Sets entrez parameters, viz. email id and API key.

- `email`: Email of the user (optional)
- `api_key`: API key of the user (optional)

This function can show the following warning

- `EmailWarning`: If no email is provided
- `ApiWarning`: If no API key is provided

**Note**

It is advised to add your email to [NCBI](https://www.ncbi.nlm.nih.gov/) before using the function. You can have your
own API key. For details about how to get API
key [click here](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us).

### `write_exome_fasta`

Creates a fasta file of exome if not exists previously or is empty

- `file_name`: The name or path of the file to be created
- `nuc_file_path`: The path of nucleotide file from where the exome is tobe extracted
- `organism_name`: Name of the organism

This function can show the following warning

- `FileNotEmptyError`: If the given file to write is not empty

### `write_nucleotide_fasta`

Creates a fasta file of nucleotides if not exists previously or is empty

- `file_name`: The name or path of the file to be created
- `cds_lst`: The tuple of FeatureLocation objects
- `record`: The SeqRecord object containing whole sequence
- `organism_name`: Name of the organism

This function can show the following warning

- `FileNotEmptyError`: If the given file to write is not empty

For details about FeatureLocation object [click here](https://biopython.org/docs/1.75/api/Bio.SeqFeature.html). For
details about SeqRecord object [click here](https://biopython.org/docs/1.75/api/Bio.SeqRecord.html).

### `write_protein_fasta`

Creates a fasta file of proteins if not exists previously or is empty

- `file_name`: The name or path of the file to be created
- `cds_lst`: The tuple of FeatureLocation objects
- `record`: The SeqRecord object containing whole sequence

This function can show the following warning

- `FileNotEmptyError`: If the given file to write is not empty

<p align="center">&copy; 2023 Souradipto Choudhuri</p>
