from os.path import join
from typing import Optional
import matplotlib.pyplot as plt
from CodonU.correspondence_analysis import build_contingency_table_aa_count, ca_aa
from CodonU.file_handler import make_dir
from CodonU.file_handler.internal_comp import is_file_writeable


def plot_ca_aa_freq_aa(handle: str, genetic_table_num: int, min_len_threshold: int = 66, n_components: int = 2,
                       organism_name: Optional[str] = None, save_image: bool = False,
                       folder_path: str = 'Report'):
    """
    Plots CA of aa frequency for genes with frequency as scale.

    :param handle: Handle to the file, or the filename as a string
    :param genetic_table_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene (optional)
    :param n_components: Components for CA (Optional)
    :param organism_name: Name of organism (optional)
    :param save_image: Options for saving the image (optional)
    :param folder_path: Folder path where image should be saved (optional)
    """
    cont_table = build_contingency_table_aa_count(handle, genetic_table_num, min_len_threshold)
    cont_table.replace(0, 0.000001, inplace=True)
    ca = ca_aa(cont_table, n_components)
    aas = ca.column_coordinates(cont_table)

    freq_vals = [sum(cont_table[aa]) for aa in aas.index.values]

    x = aas.iloc[:, 0]
    y = aas.iloc[:, 1]
    fig = plt.figure(figsize=(16, 9))
    plt.scatter(x=x, y=y, alpha=0.8, c=freq_vals, cmap='viridis', zorder=2)
    plt.grid(True, linestyle=':')
    plt.axvline(0, color='red', zorder=1)
    plt.axhline(0, color='red', zorder=1)
    plt.xlabel(f'Axis 0 (Variance: {round(ca.percentage_of_variance_[0], 4)}%)')
    plt.ylabel(f'Axis 1 (Variance: {round(ca.percentage_of_variance_[1], 4)}%)')
    for x, y, t in zip(x, y, aas.index.values):
        x = x * (1 + 0.01)
        y = y * (1 + 0.01)
        plt.text(x, y, t)
    c_bar = plt.colorbar()
    c_bar.set_label('Frequency')
    plt.title(f'Total genes: {len(cont_table.iloc[:, 0])}')
    sup_title = f'CA of AA Frequency of {organism_name} [AAs]' if organism_name else 'CA of AA Frequency [AAs]'
    plt.suptitle(sup_title)

    if save_image:
        make_dir(folder_path)
        name = f'CA_aa_freq_aa_{organism_name}.png' if organism_name else 'CA_aa_freq_aa.png'
        file_name = join(folder_path, name)
        if is_file_writeable(file_name):
            plt.savefig(file_name, dpi=500)

    plt.show()
    plt.close(fig)
