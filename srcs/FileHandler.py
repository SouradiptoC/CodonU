import os
from Bio import Entrez
from Bio.SeqIO import write
from Bio.SeqRecord import SeqRecord
from CdsExtractor import extract_cds, extract_prot
from Errors import FileNotEmpty, NoEmailError


def set_entrez_email(email: str | None) -> None:
    """
    Sets Bio.Entrez.email parameter to given email
    :param email: Email of user
    """
    if not email:
        print('Setting provided email to entrez.email')
        Entrez.email = email
    else:
        raise NoEmailError


def set_entrez_api_key(api_key: str | None) -> None:
    """
    Sets Bio.Entrez.api_key parameter to given api_key
    :param api_key: API key of the user
    """
    if not api_key:
        print('Setting provided API key to entrez.api_key')
        Entrez.api_key = api_key


def set_entrez_param(email: str | None = None, api_key: str | None = None) -> None:
    """
    Sets entrez parameters
    :param email: Email of the user
    :param api_key: API key of the user (optional)
    """
    set_entrez_email(email)
    set_entrez_api_key(api_key)


def make_dir(path: str) -> None:
    """
    Makes a directory if not present already
    :param path: Path of the directory
    """
    if not os.path.isdir(path):
        os.mkdir(path)


def is_file(path: str) -> bool:
    """
    Checks if a file exists or not
    :param path: Path of the file
    :return: True if present else False
    """
    return os.path.isfile(path)


def is_file_empty(path: str) -> bool:
    if os.stat(path).st_size == 0:
        return True
    else:
        raise FileNotEmpty(path)


def write_nucleotide_fasta(file_name: str, cds_lst: list, record: SeqRecord) -> None:
    """
    Creates a fasta file of nucleotides if not exists previously or is empty
    :param file_name: The name or path of the file
    :param cds_lst: The list of FeatureLocation objects
    :param record: The SeqRecord object containing whole sequence
    """
    if not is_file(file_name) or is_file_empty(file_name):
        with open(file_name, 'w') as out_file:
            for i in range(len(cds_lst)):
                cds = extract_cds(record, cds_lst[i], i + 1)
                write(cds, out_file, 'fasta')


def write_protein_fasta(file_name: str, cds_lst: list, organism_name: str) -> None:
    """
    Creates a fasta file of proteins if not exists previously or is empty
    :param file_name: The name or path of the file
    :param cds_lst: The list of FeatureLocation objects
    :param organism_name: Name of the organism
    """
    if not is_file(file_name) or is_file_empty(file_name):
        with open(file_name, 'w') as out_file:
            for i in range(len(cds_lst)):
                cds = extract_prot(cds_lst[i], organism_name, i + 1)
                write(cds, out_file, 'fasta')


if __name__ == '__main__':
    make_dir('Result')
