from Bio import Entrez
from Bio.SeqIO import read
from Bio.SeqRecord import SeqRecord
from CodonU.cua_logger import *
from CodonU.file_handler.network_speed import test_speed
import sys


def get_gb(accession_id: str) -> SeqRecord:
    """
    Gets the Sequence Record object from a given accession number

    :param accession_id: Provided accession number
    :return: The Sequence Record object
    """
    try:
        test_speed()
        console_log.info(f"Retrieval started")
        file_log.info(f"Retrieval started")
        handle = Entrez.efetch(db='nucleotide', id=accession_id, rettype='gb', retmode='text')
        record = read(handle, 'gb')
        console_log.info(f"Genbank file of {record.annotations['source']} retrieved successfully")
        file_log.info(f"Genbank file of {record.annotations['source']} retrieved successfully")
        return record
    except Exception as e:
        console_log.error(f'Following error occurred. See log files for details\n{e}')
        file_log.exception(e)
        sys.exit(-1)
