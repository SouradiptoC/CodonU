from typing import Any
from warnings import filterwarnings
from Bio.SeqIO import parse
from .plot_funcs import _plot_pr2
from CodonU.analyzer.internal_comp import filter_reference, gc_123, at_123, g3, a3


def plot_pr2(handle: str | Any, min_len_threshold: int, organism_name: str | None = None, save_image: bool = False,
             folder_path: str = '', gene_analysis: bool = True):
    """
    Plots A3/AT3 values against G3/GC3 values from given fasta file

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
    gc3_val_lst = []
    g3_val_lst = []
    at3_val_lst = []
    a3_val_lst = []
    for seq in reference:
        _, _, _, gc3 = gc_123(seq)
        _, _, _, at3 = at_123(seq)
        g_3 = g3(seq)
        a_3 = a3(seq)
        gc3_val_lst.append(gc3 / 100)
        at3_val_lst.append(at3 / 100)
        g3_val_lst.append(g_3 / 100)
        a3_val_lst.append(a_3 / 100)

    _plot_pr2(gc3_val_lst, at3_val_lst, g3_val_lst, a3_val_lst, organism_name, save_image, folder_path, gene_analysis)
