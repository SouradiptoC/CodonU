from Bio.SeqIO import write
from Bio.SeqRecord import SeqRecord
from .internal_comp import is_file, is_file_empty
from CodonU.extractor import extract_cds


def write_nucleotide_fasta(file_name: str, cds_lst: tuple, record: SeqRecord, organism_name: str) -> None:
    """
    Creates a fasta file of nucleotides if not exists previously or is empty

    :param file_name: The name of the file
    :param cds_lst: The tuple of FeatureLocation objects
    :param record: The SeqRecord object containing whole sequence
    :param organism_name: Name of the organism
    :raises FileNotEmptyError: If the given file to write is not empty
    """
    if not is_file(file_name) or is_file_empty(file_name):
        with open(file_name, 'w') as out_file:
            for i in range(len(cds_lst)):
                cds = extract_cds(record, cds_lst[i], i + 1)
                write(cds, out_file, 'fasta')
    print(f"Nucleotide file for {organism_name} created successfully")
