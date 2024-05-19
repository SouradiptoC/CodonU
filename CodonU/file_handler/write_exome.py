import os.path
import sys

from CodonU.file_handler.internal_comp import _write_exome
from CodonU.file_handler import make_dir
from CodonU.file_handler import get_gb
from CodonU.extractor import extract_cds, extract_cds_lst
from CodonU.cua_logger import *


def write_exome(accession_id: str, folder_path: str = 'Exome', exclude_stops: bool = True):
    """
    Creates a fasta file of all exones if not exists previously or is empty

    *Note:* It will ask if you want to re-write an existing file.

    :param accession_id: Accession id of organism
    :param folder_path: Intended folder path [Default to Exome]
    :param exclude_stops: If true, intermediate stops codons are excluded from exome
    """
    records = get_gb(accession_id)
    cds_feature_lst = extract_cds_lst(records[accession_id])
    cds_lst = [extract_cds(records[accession_id], cds_feature) for cds_feature in cds_feature_lst]
    make_dir(folder_path)
    file_path = os.path.join(folder_path, f'{records[accession_id].id}.fna')
    try:
        console_log.info(f'Started writing exome file of {accession_id}')
        file_log.info(f'Started writing exome file of {accession_id}')
        msg = _write_exome(file_path, cds_lst, exclude_stops, multi=False)
        console_log.info(msg)
        file_log.info(msg)
    except Exception as e:
        console_log.error(f'Following error occurred. See log files for details\n{e}')
        file_log.exception(e)
        sys.exit(-1)
