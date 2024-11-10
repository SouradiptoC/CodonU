import sys
from typing import Any

from Bio.SeqRecord import SeqRecord

from CodonU.cua_logger import *


def extract_cds_lst(record: SeqRecord) -> tuple[Any, ...]:
    """
    Extracts the list of features if their type is CDS

    :param record: Original Sequence Record object from where the CDS is to be extracted
    :return: A tuple of FeatureLocation objects
    """
    try:
        cds_lst = [cds for cds in record.features if cds.type == 'CDS' and 'pseudo' not in cds.qualifiers.keys()]
        return tuple(cds_lst)
    except Exception as e:
        console_log.error(f'Following exception occurred.\n{e}')
        file_log.exception(f'Exception occurred.\n{e}')
        sys.exit(-1)
