from Bio import Entrez
from Bio.SeqIO import read
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys

from CodonU.cua_logger import *

from CodonU.file_handler.network_speed import test_speed


def _get_gb(accession_id):
    """
    Gets the Sequence Record object from a given accession number

    :param accession_id: Provided accession number
    :return: The Sequence Record object
    """
    handle = Entrez.efetch(db='nucleotide', id=accession_id, rettype='gb', retmode='text')
    record = read(handle, 'gb')
    return record


def get_gbs(accession_ids: list[str]) -> dict[str, SeqRecord]:
    """
    Gets the Sequence Record object from a given accession number

    :param accession_ids: Provided multiple accession ids
    :return: The dict of Sequence Record object
    """
    try:
        test_speed()
        console_log.info(f"Retrieval started")
        file_log.info(f"Retrieval started")
        records = dict()
        with ThreadPoolExecutor(thread_name_prefix='codonu_genbank_retr') as executor:
            future_records = {executor.submit(_get_gb, accession_id): accession_id for accession_id in accession_ids}
            for future_record in as_completed(future_records):
                acs_id = future_records[future_record]
                try:
                    records.update({acs_id: future_record.result()})
                except Exception as e:
                    console_log.error(f"Can't fetch {acs_id}. Raising RuntimeError.")
                    file_log.exception(
                        f"Can't fetch {acs_id}. Raising RuntimeError as side effect. Initial error is\n{e}")
                    raise RuntimeError('Bad accession id')
        console_log.info(f"Genbank file of {accession_ids} retrieved successfully")
        file_log.info(f"Genbank file of {accession_ids} retrieved successfully")
        return records
    except Exception as e:
        console_log.error(f'Following error occurred. See log files for details\n{e}')
        file_log.exception(e)
        sys.exit(-1)
