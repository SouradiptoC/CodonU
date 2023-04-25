from Bio.SeqIO import read, write
from os.path import join, abspath
from CodonU.file_handler import make_dir
from CodonU.file_handler.internal_comp import is_file_writeable


def generate_phylo_input(file_lst: list[str], file_name: str = 'phylo_input', folder_path: str = 'Report'):
    """
    Generates the alignment input file

    :param file_lst: List containing paths of fasta file
    :param file_name: Output filename
    :param folder_path: Folder name
    """
    make_dir(folder_path)
    file_path = join(folder_path, file_name) + '.fasta'
    if is_file_writeable(file_path):
        for in_file in file_lst:
            records = read(in_file, 'fasta')
            write(records, open(file_path, 'a'), 'fasta')
        print(f'The alignment input file can be found at {abspath(file_path)}')
