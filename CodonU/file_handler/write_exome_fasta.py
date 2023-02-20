from Bio.SeqIO import write
from .internal_comp import is_file, is_file_empty
from CodonU.extractor import extract_exome


def write_exome_fasta(file_name: str, nuc_file_path: str, organism_name: str) -> None:
    """
    Creates a fasta file of exome if not exists previously or is empty

    :param file_name: The name of the file
    :param nuc_file_path: The path of nucleotide file
    :param organism_name: Name of the organism
    :raises FileNotEmptyError: If the given file to write is not empty
    """
    if not is_file(file_name) or is_file_empty(file_name):
        with open(file_name, 'w') as out_file:
            exome = extract_exome(nuc_file_path, organism_name)
            write(exome, out_file, 'fasta')
    print(f"Exome file for {organism_name} created successfully")
