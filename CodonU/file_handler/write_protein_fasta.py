from os.path import abspath
from Bio.SeqIO import write
from CodonU.file_handler import get_gb
from CodonU.file_handler.internal_comp import is_file, is_file_empty
from CodonU.extractor import extract_prot, extract_cds_lst


def write_protein_fasta(accession_id: str, file_path: str):
    """
    Creates a fasta file of proteins if not exists previously or is empty

    :param accession_id: Accession id of organism
    :param file_path: Intended file path
    """
    records = get_gb(accession_id)
    cds_feature_lst = extract_cds_lst(records)

    if not is_file(file_path) or is_file_empty(file_path):
        with open(file_path, 'w') as out_file:
            prot_lst = [extract_prot(cds_feature, records.annotations['organism']) for cds_feature in cds_feature_lst]
            write(prot_lst, out_file, 'fasta')
            print(f"Protein file can be found at: {abspath(file_path)}")
