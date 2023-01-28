import os
import pandas as pd
from Bio import Entrez
from Bio.SeqIO import write, read
from Bio.SeqRecord import SeqRecord
from Extractor import extract_cds, extract_prot, extract_exome
from Errors import FileNotEmptyError, EmailWarning, ApiWarning


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


def write_nucleotide_fasta(file_name: str, cds_lst: tuple, record: SeqRecord, organism_name: str) -> None:
    """
    Creates a fasta file of nucleotides if not exists previously or is empty
    :param file_name: The name of the file
    :param cds_lst: The tuple of FeatureLocation objects
    :param record: The SeqRecord object containing whole sequence
    :param organism_name: Name of the organism
    """
    if not is_file(file_name) or is_file_empty(file_name):
        with open(file_name, 'w') as out_file:
            for i in range(len(cds_lst)):
                cds = extract_cds(record, cds_lst[i], i + 1)
                write(cds, out_file, 'fasta')
    print(f"Nucleotide file for {organism_name} created successfully")


def write_protein_fasta(file_name: str, cds_lst: tuple, organism_name: str) -> None:
    """
    Creates a fasta file of proteins if not exists previously or is empty
    :param file_name: The name of the file
    :param cds_lst: The tuple of FeatureLocation objects
    :param organism_name: Name of the organism
    """
    if not is_file(file_name) or is_file_empty(file_name):
        with open(file_name, 'w') as out_file:
            for i in range(len(cds_lst)):
                cds = extract_prot(cds_lst[i], organism_name, i + 1)
                write(cds, out_file, 'fasta')
    print(f"Protein file for {organism_name} created successfully")


def write_exome_fasta(file_name: str, nuc_file_path: str, organism_name: str):
    """
    Creates a fasta file of exome if not exists previously or is empty
    :param file_name: The name of the file
    :param nuc_file_path: The path of nucleotide file
    :param organism_name: Name of the organism
    """
    if not is_file(file_name) or is_file_empty(file_name):
        with open(file_name, 'w') as out_file:
            # TODO provide threshold value
            exome = extract_exome(nuc_file_path, organism_name)
            write(exome, out_file, 'fasta')
    print(f"Exome file for {organism_name} created successfully")


if __name__ == '__main__':
    # print(is_file_empty('temp2.py'))
    write_exome_fasta('../Results/Exons/temp_2.fasta', '../Results/Nucleotide/Staphylococcus_agnetis_nucleotide.fasta',
                      'Staphylococcus agnetis')
