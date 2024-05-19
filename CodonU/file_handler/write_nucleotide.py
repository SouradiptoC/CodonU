import os
import sys

from CodonU.cua_logger import *
from CodonU.file_handler import get_gb
from CodonU.file_handler.internal_comp import _write_nucleotide


def write_nucleotide(accession_id: str, folder_path: str = 'Nuc'):
    """
    Creates a fasta file of nucleotides if not exists previously or is empty

    :param accession_id: Accession id of organism
    :param folder_path: Intended folder path (default to 'Nuc')
    """
    try:
        console_log.info(f'Started writing nucleotide file of {accession_id}')
        file_log.info(f'Started writing nucleotide file of {accession_id}')
        records = get_gb(accession_id)
        file_path = os.path.join(folder_path, f'{records[accession_id].id}.ffn')
        msg = _write_nucleotide(file_path, records[accession_id], multi=False)
        console_log.info(msg)
        file_log.info(msg)
    except Exception as e:
        console_log.error(f'Following error occurred. See log files for details\n{e}')
        file_log.exception(e)
        sys.exit(-1)
