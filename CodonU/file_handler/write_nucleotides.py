import os
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from CodonU.cua_logger import *
from CodonU.file_handler import get_gbs, make_dir
from CodonU.file_handler.internal_comp import _write_nucleotide


def write_nucleotides(accession_ids: list[str], folder_path: str = 'Nuc'):
    """
    Creates a fasta file of nucleotides if not exists previously or is empty

    :param accession_ids: List of accession ids of organisms
    :param folder_path: Intended folder path (default to 'Nuc')
    """
    records = get_gbs(accession_ids)
    file_path_dict = {accession_id: os.path.join(folder_path, f'{records[accession_id].id}.ffn') for accession_id in
                      accession_ids}
    make_dir(folder_path)
    console_log.info(f'Started writing nucleotide file of {accession_ids}')
    file_log.info(f'Started writing nucleotide file of {accession_ids}')
    with ThreadPoolExecutor(thread_name_prefix='codonu_write_nucleo') as executor:
        future_records = {
            executor.submit(_write_nucleotide, file_path_dict[acs_id], records[acs_id], True): acs_id
            for acs_id in accession_ids
        }
        for future_record in as_completed(future_records):
            acs_id = future_records[future_record]
            try:
                msg = future_record.result()
                console_log.info(msg)
                file_log.info(msg)
            except Exception as e:
                console_log.error(f'Following exception occurred during writing file of {acs_id}.\n{e}')
                file_log.exception(f'Exception occurred during writing file of {acs_id}.\n{e}')
                sys.exit(-1)
