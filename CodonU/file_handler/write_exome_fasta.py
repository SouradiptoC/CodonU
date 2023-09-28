from os.path import abspath
from Bio.SeqIO import write
from CodonU.file_handler.internal_comp import is_file, is_file_empty
from CodonU.extractor import extract_exome
from CodonU.file_handler import get_gb
from CodonU.extractor import extract_cds, extract_cds_lst


def write_exome_fasta(accession_id: str, file_path: str, exclude_stops: bool = True):
    """
    Creates a fasta file of all exones if not exists previously or is empty

    :param accession_id: Accession id of organism
    :param file_path: Intended file path
    :param exclude_stops: If true, intermediate stops codons are excluded from exome
    :return:
    """
    records = get_gb(accession_id)
    cds_feature_lst = extract_cds_lst(records)
    cds_lst = [extract_cds(records, cds_feature) for cds_feature in cds_feature_lst]
    if not is_file(file_path) or is_file_empty(file_path):
        with open(file_path, 'w') as out_file:
            exome = extract_exome(cds_lst, exclude_stops)
            write(exome, out_file, 'fasta')
        print(f"Exome file can be found at: {abspath(file_path)}")
