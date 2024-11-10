import sys

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from CodonU.extractor.internal_comp import extract_cds_seq
from CodonU.cua_logger import *


def extract_cds(record: SeqRecord, cds_feature: SeqFeature) -> SeqRecord:
    """
    Returns the CDS as a Sequence Record object

    :param record: Original Sequence Record object from where the CDS is to be extracted
    :param cds_feature: Sequence Feature object of corresponding SeqRecord object
    :return: The new Sequence Record object containing the CDS
    """
    try:
        desc_lst = [f'[{key}={val}]' for key, val in cds_feature.qualifiers.items() if not key == 'translation']
        cds = SeqRecord(
            seq=extract_cds_seq(record.seq, cds_feature.location),
            id=f"{record.id}|{record.annotations['organism']}",
            name=f"{record.annotations['organism']}",
            description=' '.join(desc_lst),
        )
        return cds
    except Exception as e:
        console_log.error(f'Following exception occurred.\n{e}')
        file_log.exception(f'Exception occurred.\n{e}')
        sys.exit(-1)
