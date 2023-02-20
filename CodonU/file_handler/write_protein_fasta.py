from Bio.SeqIO import write
from .internal_comp import is_file, is_file_empty
from CodonU.extractor import extract_prot


def write_protein_fasta(file_name: str, cds_lst: tuple, organism_name: str) -> None:
    """
    Creates a fasta file of proteins if not exists previously or is empty

    :param file_name: The name of the file
    :param cds_lst: The tuple of FeatureLocation objects
    :param organism_name: Name of the organism
    :raises FileNotEmptyError: If the given file to write is not empty
    """
    if not is_file(file_name) or is_file_empty(file_name):
        with open(file_name, 'w') as out_file:
            for i in range(len(cds_lst)):
                cds = extract_prot(cds_lst[i], organism_name, i + 1)
                write(cds, out_file, 'fasta')
    print(f"Protein file for {organism_name} created successfully")
