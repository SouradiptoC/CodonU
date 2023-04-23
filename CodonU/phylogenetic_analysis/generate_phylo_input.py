from Bio.SeqIO import read, write
from os.path import join, abspath, isfile
from os import stat
from CodonU.cua_errors import FileNotEmptyError
from CodonU.file_handler import make_dir


def _is_file_writeable(path: str) -> bool:
    """
    Checks if provided file exists and is writable

    :param path: Path to the file
    :raise FileNotEmptyError: If user doesn't want to append the fasta sequence
    :return: True if user wants to, False otherwise
    """
    if isfile(path) and stat(path).st_size != 0:
        flg = input('Provided file not empty! Sequences will be appended to the file. Proceed [y/n]?: ')
        if flg in ['y', 'Y']:
            return True
        raise FileNotEmptyError(path)
    else:
        return True


def generate_phylo_input(file_lst: list[str], file_name: str = 'phylo_input', folder_path: str = 'Report'):
    """
    Generates the alignment input file

    :param file_lst: List containing paths of fasta file
    :param file_name: Output filename
    :param folder_path: Folder name
    """
    make_dir(folder_path)
    file_path = join(folder_path, file_name) + '.fasta'
    if _is_file_writeable(file_path):
        for in_file in file_lst:
            records = read(in_file, 'fasta')
            write(records, open(file_path, 'a'), 'fasta')
        print(f'The alignment input file can be found at {abspath(file_path)}')
