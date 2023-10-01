from os.path import join
from typing import Optional
import matplotlib.pyplot as plt
from Bio.SeqIO import parse
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from CodonU.correspondence_analysis import build_contingency_table_aa_count, ca_aa
from CodonU.file_handler import make_dir
from CodonU.file_handler.internal_comp import is_file_writeable
from CodonU.analyzer import filter_reference


def plot_ca_aa_freq_gene(handle: str, genetic_table_num: int, scale: str, min_len_threshold: int = 66,
                         n_components: int = 2, organism_name: Optional[str] = None, save_image: bool = False,
                         folder_path: str = 'Report'):
    """
    Plots CA for genes with aromaticity or gravy score as scale.\n
    **NOTE** Values of `scale` supported are `'gravy'`, `'aroma'`

    :param handle: Handle to the file, or the filename as a string
    :param genetic_table_num: Genetic table number for codon table
    :param scale: Scale for gene colors
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene (optional)
    :param n_components: Components for CA (Optional)
    :param organism_name: Name of organism (optional)
    :param save_image: Options for saving the image (optional)
    :param folder_path: Folder path where image should be saved (optional)
    """
    cont_table = build_contingency_table_aa_count(handle, genetic_table_num, min_len_threshold)
    cont_table.replace(0, 0.000001, inplace=True)
    ca = ca_aa(cont_table, n_components)
    genes = ca.row_coordinates(cont_table)

    if scale == 'aroma':
        color_bar = [ProteinAnalysis(str(prot.seq)).aromaticity() for prot in
                     filter_reference(parse(handle, 'fasta'), min_len_threshold, 'aa')]
        label = 'Aroma Values'
    elif scale == 'gravy':
        color_bar = [ProteinAnalysis(str(prot.seq)).gravy() for prot in
                     filter_reference(parse(handle, 'fasta'), min_len_threshold, 'aa')]
        label = 'GRAVY values'
    else:
        raise ValueError(f"Given scale type {scale} not supported. Supported scales are 'gravy' and 'aroma'.")

    x = genes.iloc[:, 0]
    y = genes.iloc[:, 1]
    fig = plt.figure(figsize=(16, 9))
    plt.scatter(x=x, y=y, alpha=0.8, c=color_bar, cmap='viridis', zorder=2)
    plt.grid(True, linestyle=':')
    plt.axvline(0, color='red', zorder=1)
    plt.axhline(0, color='red', zorder=1)
    plt.xlabel(f'Axis 0 (Variance: {round(ca.percentage_of_variance_[0], 4)}%)')
    plt.ylabel(f'Axis 1 (Variance: {round(ca.percentage_of_variance_[1], 4)}%)')
    c_bar = plt.colorbar()
    c_bar.set_label(label)
    plt.title(f'Total genes: {len(cont_table.iloc[:, 0])}')
    sup_title = f'CA of AA Frequency of {organism_name} [genes]' if organism_name else 'CA of AA Frequency [genes]'
    plt.suptitle(sup_title)

    if save_image:
        make_dir(folder_path)
        name = f'CA_aa_freq_aa_{organism_name}.png' if organism_name else 'CA_aa_freq_aa.png'
        file_name = join(folder_path, name)
        if is_file_writeable(file_name):
            plt.savefig(file_name, dpi=500)

    plt.show()
    plt.close(fig)
