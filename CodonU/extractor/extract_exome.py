import sys

from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord

from CodonU.cua_logger import *


def extract_exome(records: list[SeqRecord], exclude_stops: bool = True) -> SeqRecord:
    """
    Extracts the exome from given nucleotides
    :param records: List of SeqRecord objects containing each CDS
    :param exclude_stops: If true, intermediate stops codons are excluded from exome
    :return: The exome
    """
    try:
        m_seq = MutableSeq('')
        tot_len = len(records)
        if exclude_stops:
            for idx, record in enumerate(records, start=1):
                m_seq += record.seq[:-3]
                if idx == tot_len:
                    m_seq += record.seq[-3:]
        else:
            for record in records:
                m_seq += record.seq

        exome = SeqRecord(
            seq=m_seq,
            id=f"{records[0].id}",
            name=f"{records[0].id.split('|')[-1]}",
            description=f"{records[0].id.split('|')[-1]} complete genome"
        )
        return exome
    except Exception as e:
        console_log.error(f'Following exception occurred.\n{e}')
        file_log.exception(f'Exception occurred.\n{e}')
        sys.exit(-1)
