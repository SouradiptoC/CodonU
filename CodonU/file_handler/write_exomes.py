import os.path
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

from CodonU.file_handler.internal_comp import _write_exome
from CodonU.file_handler import make_dir
from CodonU.file_handler import get_gbs
from CodonU.extractor import extract_cds, extract_cds_lst
from CodonU.cua_logger import *


def write_exomes(accession_ids: list[str], folder_path: str = 'Exome', exclude_stops: bool = True):
    """
    Creates fasta files of all exones of all organisms

    *Note:* It will ask if you want to re-write an existing file.

    :param accession_ids: List of accession id of organism
    :param folder_path: Intended folder path [Default to Exome]
    :param exclude_stops: If true, intermediate stops codons are excluded from exome
    """
    records = get_gbs(accession_ids)
    cds_feature_dict = {accession_id: extract_cds_lst(records[accession_id]) for accession_id in accession_ids}
    cds_dic = {accession_id: [extract_cds(records[accession_id], cds_feature_lst) for cds_feature_lst in
                              cds_feature_dict[accession_id]] for accession_id in accession_ids}
    make_dir(folder_path)
    file_path_dict = {accession_id: os.path.join(folder_path, f'{records[accession_id].id}.fna') for accession_id in
                      accession_ids}
    console_log.info(f'Started writing exome file of {accession_ids}')
    file_log.info(f'Started writing exome file of {accession_ids}')
    with ThreadPoolExecutor(thread_name_prefix='codonu_write_exome') as executor:
        future_records = {
            executor.submit(_write_exome, file_path_dict[acs_id], cds_dic[acs_id], exclude_stops, True): acs_id
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
