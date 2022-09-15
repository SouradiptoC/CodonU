import os
from Bio.SeqIO import write
from Bio.SeqRecord import SeqRecord
from CdsExtractor import extract_cds


def make_dir(path: str) -> None:
    """
    Makes a directory if not present already
    :param path: Path of the directory
    """
    if not os.path.isdir(path):
        os.mkdir(path)


def is_file(path: str) -> bool:
    """
    Checks if a file exists or not
    :param path: Path of the file
    :return: True if present else False
    """
    return os.path.isfile(path)


def write_fasta(file_name: str, cds_lst: list, record: SeqRecord) -> None:
    """
    Creates a fasta file if not exists previously
    :param file_name: The name or path of the file
    :param cds_lst: The list of FeatureLocation object
    :param record: The SeqRecord object containing whole sequence
    """
    if not is_file(file_name):
        with open(file_name, 'a') as out_file:
            for i in range(len(cds_lst)):
                cds = extract_cds(record, cds_lst[i], i + 1)
                write(cds, out_file, 'fasta')


if __name__ == '__main__':
    make_dir('Result')
