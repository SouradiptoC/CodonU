from os.path import join
from typing import Optional
import matplotlib.pyplot as plt
from CodonU.correspondence_analysis import build_contingency_table_codon_count, ca_codon
from CodonU.file_handler import make_dir
from CodonU.file_handler.internal_comp import is_file_writeable


def plot_ca_codon_freq_codon(handle: str, genetic_table_num: int, min_len_threshold: int = 200, n_components: int = 2,
                             single_syn_codons: Optional[list] = None, organism_name: Optional[str] = None,
                             save_image: bool = False, folder_path: str = 'Report'):
    """
    Plots CA of codon frequency for codons with frequency as scale

    :param handle: Handle to the file, or the filename as a string
    :param genetic_table_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene (optional)
    :param n_components: Components for CA (Optional)
    :param single_syn_codons: List of codons belonging to SF1 (optional)
    :param organism_name: Name of organism (optional)
    :param save_image: Options for saving the image (optional)
    :param folder_path: Folder path where image should be saved (optional)
    """
    if single_syn_codons is None:
        single_syn_codons = ['ATG', 'TGG']

    cont_table = build_contingency_table_codon_count(handle, genetic_table_num, min_len_threshold)
    cont_table.replace(0, 0.000001, inplace=True)
    ca = ca_codon(cont_table, n_components, single_syn_codons)
    cont_table.drop(single_syn_codons, axis=1, inplace=True)
    codons = ca.column_coordinates(cont_table)

    freq_vals = [sum(cont_table[codon]) for codon in codons.index.values]

    x = codons.iloc[:, 0]
    y = codons.iloc[:, 1]
    fig = plt.figure(figsize=(16, 9))
    plt.scatter(x=x, y=y, alpha=0.8, c=freq_vals, cmap='viridis', zorder=2)
    plt.grid(True, linestyle=':')
    plt.axvline(0, color='red', zorder=1)
    plt.axhline(0, color='red', zorder=1)
    plt.xlabel(f'Axis 0 (Variance: {round(ca.percentage_of_variance_[0], 4)}%)')
    plt.ylabel(f'Axis 1 (Variance: {round(ca.percentage_of_variance_[1], 4)}%)')
    for x, y, t in zip(x, y, codons.index.values):
        x = x * (1 + 0.01)
        y = y * (1 + 0.01)
        plt.text(x, y, t)
    c_bar = plt.colorbar()
    c_bar.set_label('Frequency')
    plt.title(f'Total genes: {len(cont_table.iloc[:, 0])}')
    sup_title = f'CA of Codon Frequency of {organism_name} [codons]' if organism_name \
        else 'CA of Codon Frequency [codons]'
    plt.suptitle(sup_title)

    if save_image:
        make_dir(folder_path)
        name = f'CA_codon_freq_codon_{organism_name}.png' if organism_name else 'CA_codon_freq_codon.png'
        file_name = join(folder_path, name)
        if is_file_writeable(file_name):
            plt.savefig(file_name, dpi=500)

    plt.show()
    plt.close(fig)
