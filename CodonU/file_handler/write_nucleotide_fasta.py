from Bio.SeqIO import write
from CodonU.file_handler.internal_comp import is_file, is_file_empty
from CodonU.extractor import extract_cds, extract_cds_lst
from CodonU.file_handler import get_gb
from os.path import abspath


def write_nucleotide_fasta(accession_id: str, file_path: str):
    """
    Creates a fasta file of nucleotides if not exists previously or is empty

    :param accession_id: Accession id of organism
    :param file_path: Intended file path
    """
    records = get_gb(accession_id)
    cds_feature_lst = extract_cds_lst(records)

    if not is_file(file_path) or is_file_empty(file_path):
        with open(file_path, 'w') as out_file:
            cds_lst = [extract_cds(records, cds_feature) for cds_feature in cds_feature_lst]
            write(cds_lst, out_file, 'fasta')
            print(f"Nucleotide file can be found at: {abspath(file_path)}")
