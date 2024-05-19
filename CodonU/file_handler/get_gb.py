from Bio.SeqRecord import SeqRecord
import sys

from CodonU.cua_logger import *
from CodonU.file_handler.network_speed import test_speed
from CodonU.file_handler.internal_comp import _get_gb


def get_gb(accession_id: str) -> dict[str, SeqRecord]:
    """
    Gets the Sequence Record object from a given accession number

    :param accession_id: Provided accession number
    :return: Dictionary with accession id as key, SeqRecord as val
    """
    try:
        test_speed()
        console_log.info(f"Retrieval started")
        file_log.info(f"Retrieval started")
        record = {accession_id: _get_gb(accession_id)}
        console_log.info(f"Genbank file of accession id: {accession_id} retrieved successfully")
        file_log.info(f"Genbank file of accession id: {accession_id} retrieved successfully")
        return record
    except Exception as e:
        console_log.error(f'Following error occurred. See log files for details\n{e}')
        file_log.exception(e)
        sys.exit(-1)
