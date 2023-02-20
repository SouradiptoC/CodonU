from typing import Any
from warnings import filterwarnings
from Bio.SeqIO import parse
from .plot_funcs import _plot_neutrality
from CodonU.analyzer.internal_comp import gc_123, filter_reference


def plot_neutrality(handle: str | Any, min_len_threshold: int, organism_name: str | None = None,
                    save_image: bool = False, folder_path: str = '', gene_analysis: bool = True):
    """
    Plots neutrality plot from given fasta file

    :param handle: Handle to the file, or the filename as a string
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :param organism_name: Name of organism (optional)
    :param save_image: Options for saving the image (optional)
    :param folder_path: Folder path where image should be saved (optional)
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    """
    filterwarnings('ignore')
    records = parse(handle, 'fasta')
    reference = filter_reference(records, min_len_threshold)
    gc_12_lst = []
    gc_3_lst = []
    for seq in reference:
        _, gc_1, gc_2, gc_3 = gc_123(seq)
        # taking avg of gc_1 and gc_2
        gc_12_lst.append((gc_1 + gc_2) / 2 / 100)
        gc_3_lst.append(gc_3 / 100)

    _plot_neutrality(gc_12_lst, gc_3_lst, organism_name, save_image, folder_path, gene_analysis)
