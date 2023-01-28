import os
import pandas as pd
import warnings
from Bio import Entrez
from Bio.SeqIO import write, read, parse
from Bio.SeqRecord import SeqRecord
from Extractor import extract_cds, extract_prot, extract_exome
from Errors import FileNotEmptyError, NoEmailError, EmailWarning, ApiWarning


def set_entrez_email(email: str | None) -> None:
    """
    Sets Bio.Entrez.email parameter to given email
    :param email: Email of user
    """
    if email:
        print('Setting provided email to entrez.email')
        Entrez.email = email
    else:
        warning = EmailWarning()
        warning.warn()


def set_entrez_api_key(api_key: str | None) -> None:
    """
    Sets Bio.Entrez.api_key parameter to given api_key
    :param api_key: API key of the user
    """
    if api_key:
        print('Setting provided API key to entrez.api_key')
        Entrez.api_key = api_key
    else:
        warning = ApiWarning()
        warning.warn()


def set_entrez_param(email: str | None = None, api_key: str | None = None) -> None:
    """
    Sets entrez parameters
    :param email: Email of the user
    :param api_key: API key of the user (optional)
    """
    set_entrez_email(email)
    set_entrez_api_key(api_key)


def get_gb(accession_id: str) -> SeqRecord:
    """
    Gets the Sequence Record object from a given accession number
    :param accession_id: Provided accession number
    :return: The Sequence Record object
    """
    print(f"Retrieval started")
    handle = Entrez.efetch(db='nucleotide', id=accession_id, rettype='gb', retmode='text')
    record = read(handle, 'gb')
    print(f"Genbank file of {record.annotations['source']} retrieved successfully")
    return record


def make_dir(path: str) -> None:
    """
    Makes a directory if not present already
    :param path: Path of the directory
    """
    if not os.path.isdir(path):
        os.mkdir(path)
        name = path.split('/')
        print(f"{name[-1]} created successfully")


def is_file_empty(path: str) -> bool:
    """
    Checks if an existing file is empty
    :param path: Path to the file
    :return: True if empty else false
    """
    if os.stat(path).st_size == 0:
        return True
    else:
        name = path.split('/')
        dec = input(f"{name[-1]} already exists. Want to re-write (y/n): ")
        if dec == 'y':
            return True
        raise FileNotEmptyError(path)


def is_file(path: str):
    return os.path.isfile(path)


def read_file(file_name: str) -> pd.DataFrame:
    """
    Returns a dataframe from given csv file
    :param file_name: Name or path to csv file
    :return: The data frame object
    """
    df = pd.read_csv(file_name)
    return df


def write_nucleotide_fasta(file_name: str, cds_lst: tuple, record: SeqRecord) -> None:
    """
    Creates a fasta file of nucleotides if not exists previously or is empty
    :param file_name: The name of the file
    :param cds_lst: The tuple of FeatureLocation objects
    :param record: The SeqRecord object containing whole sequence
    """
    if not is_file(file_name) or is_file_empty(file_name):
        with open(file_name, 'w') as out_file:
            for i in range(len(cds_lst)):
                cds = extract_cds(record, cds_lst[i], i + 1)
                write(cds, out_file, 'fasta')
    print(f"Nucleotide file for {file_name.split('/')[-1]} created successfully")


def write_protein_fasta(file_name: str, cds_lst: tuple, organism_name: str) -> None:
    """
    Creates a fasta file of proteins if not exists previously or is empty
    :param file_name: The name of the file
    :param cds_lst: The tuple of FeatureLocation objects
    :param organism_name: Name of the organism
    :return:
    """
    if not is_file(file_name) or is_file_empty(file_name):
        with open(file_name, 'w') as out_file:
            for i in range(len(cds_lst)):
                cds = extract_prot(cds_lst[i], organism_name, i + 1)
                write(cds, out_file, 'fasta')
    print(f"Protein file for {file_name.split('/')[-1]} created successfully")


def write_exome_fasta(file_name: str, nuc_file_path: str, organism_name: str):
    if not is_file(file_name) or is_file_empty(file_name):
        with open(file_name, 'w') as out_file:
            exome = extract_exome(nuc_file_path, organism_name)


if __name__ == '__main__':
    # print(is_file_empty('temp2.py'))
    x = None
    print(set_entrez_email(x))
