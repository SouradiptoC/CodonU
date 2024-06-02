import os
from typing import Optional
from Bio import Entrez
from Bio.SeqIO import read, write
from Bio.SeqRecord import SeqRecord

from CodonU.cua_errors import FileNotEmptyError
from CodonU.cua_warnings import EmailWarning, ApiWarning
from CodonU.extractor import extract_exome, extract_cds_lst, extract_cds, extract_prot
from CodonU.cua_logger import *


def set_entrez_email(email: Optional[str]) -> None:
    """
    Sets Bio.Entrez.email parameter to given email

    :param email: Email of user
    :raises EmailWarning: If no email is provided
    """
    if email:
        print('Setting provided email to entrez.email')
        Entrez.email = email
    else:
        warning = EmailWarning()
        warning.warn()


def set_entrez_api_key(api_key: Optional[str]) -> None:
    """
    Sets Bio.Entrez.api_key parameter to given api_key

    :param api_key: API key of the user
    :raises ApiWarning: If no API key is provided
    """
    if api_key:
        print('Setting provided API key to entrez.api_key')
        Entrez.api_key = api_key
    else:
        warning = ApiWarning()
        warning.warn()


def is_file_empty(path: str) -> bool:
    """
    Checks if an existing file is empty

    :param path: Path to the file
    :return: True if empty else false
    :raises FileNotEmptyError: If the given file to write is not empty
    """
    if os.stat(path).st_size == 0:
        return True
    else:
        name = path.split('/')
        dec = input(f"{name[-1]} already exists. Want to re-write (y/n): ")
        if dec == 'y':
            return True
        raise FileNotEmptyError(path)


def is_file(path: str) -> bool:
    """
    Checks if file exists or not

    :param path: Path to the file
    :return: True if exists else False
    """
    return os.path.isfile(path)


def is_file_writeable(path: str):
    if is_file(path) and os.stat(path).st_size != 0:
        flg = input(
            'Provided file not empty! Your action will result into completely changing the content of the file. '
            'Proceed [y/n]?: ')
        if flg in ['y', 'Y']:
            return True
        raise FileNotEmptyError(path)
    else:
        return True


def _get_gb(accession_id):
    """
    Gets the Sequence Record object from a given accession number

    :param accession_id: Provided accession number
    :return: The Sequence Record object
    """
    # warnings.filterwarnings('ignore')
    handle = Entrez.efetch(db='nucleotide', id=accession_id, rettype='gb', retmode='text')
    record = read(handle, 'gb')
    # warnings.resetwarnings()
    return record


def _write_exome(file_path: str, cds_lst: list[SeqRecord], exclude_stops: bool = True, multi: bool = False) -> str:
    """
    Creates a fasta file of all exones (if multi is false then asks if user want to re-write)

    :param file_path: File path to write at
    :param cds_lst: CDS list of organism
    :param exclude_stops: If true, intermediate stops codons are excluded from exome
    :param multi: If true, means the function is called from multi-threading
    :return: File path of created file
    :raises RuntimeError: In case of any exception
    """
    try:
        if multi or (not is_file(file_path) or is_file_empty(file_path)):
            with open(file_path, 'w') as out_file:
                exome = extract_exome(cds_lst, exclude_stops)
                write(exome, out_file, 'fasta')
            return f'Exome file can be found at: {os.path.abspath(file_path)}'
    except FileNotEmptyError as fne:
        console_log.error(f'Following error occurred. See log files for details\n{fne}')
        file_log.exception(fne)
        raise RuntimeError('Bad file')
    except Exception as exe:
        console_log.error(f'Following error occurred. See log files for details\n{exe}')
        file_log.exception(exe)
        raise RuntimeError


def _write_nucleotide(file_path: str, record: SeqRecord, multi: bool = False):
    """
    Creates .ffn file for genetic code
    :param file_path: file path to write
    :param record: SeqRecord obj
    :param multi: If true, being called from multi-threading
    :return: path to file
    :raises RuntimeError: If file not empty or else
    """
    try:
        cds_feature_lst = extract_cds_lst(record)
        if multi or (not is_file(file_path) or is_file_empty(file_path)):
            with open(file_path, 'w') as out_file:
                cds_lst = [extract_cds(record, cds_feature) for cds_feature in cds_feature_lst]
                write(cds_lst, out_file, 'fasta')
                return f'Nucleotide file can be found at: {os.path.abspath(file_path)}'
    except FileNotEmptyError as fne:
        console_log.error(f'Following error occurred. See log files for details\n{fne}')
        file_log.exception(fne)
        raise RuntimeError('Bad file')
    except Exception as exe:
        console_log.error(f'Following error occurred. See log files for details\n{exe}')
        file_log.exception(exe)
        raise RuntimeError


def _write_protein(file_path: str, record: SeqRecord, multi: bool = False):
    """
    Creates .faa file for genetic code
    :param file_path: file path to write
    :param record: SeqRecord obj
    :param multi: If true, being called from multi-threading
    :return: path to file
    :raises RuntimeError: If file not empty or else
    """
    try:
        cds_feature_lst = extract_cds_lst(record)

        if multi or (not is_file(file_path) or is_file_empty(file_path)):
            with open(file_path, 'w') as out_file:
                prot_lst = [extract_prot(cds_feature, record.annotations['organism']) for cds_feature in
                            cds_feature_lst]
                write(prot_lst, out_file, 'fasta')
                return f'Protein file can be found at: {os.path.abspath(file_path)}'
    except FileNotEmptyError as fne:
        console_log.error(f'Following error occurred. See log files for details\n{fne}')
        file_log.exception(fne)
        raise RuntimeError('Bad file')
    except Exception as exe:
        console_log.error(f'Following error occurred. See log files for details\n{exe}')
        file_log.exception(exe)
        raise RuntimeError
