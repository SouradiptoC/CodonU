from os.path import join
from typing import Optional
import matplotlib.pyplot as plt
from Bio.SeqIO import parse
from CodonU.analyzer import filter_reference
from CodonU.correspondence_analysis import build_contingency_table_codon_rscu, ca_codon
from CodonU.file_handler import make_dir
from CodonU.file_handler.internal_comp import is_file_writeable


def plot_ca_codon_rscu_gene(handle: str, genetic_table_num: int, min_len_threshold: int = 200, n_components: int = 2,
                            single_syn_codons: Optional[list] = None, organism_name: Optional[str] = None,
                            save_image: bool = False, folder_path: str = 'Report'):
    """
    Plots CA of codon RSCU for codons with frequency as scale

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

    cont_table = build_contingency_table_codon_rscu(handle, genetic_table_num, min_len_threshold)
    cont_table.replace(0, 0.000001, inplace=True)
    ca = ca_codon(cont_table, n_components, single_syn_codons)
    cont_table.drop(single_syn_codons, axis=1, inplace=True)
    # converting to mitigate TypeError: no supported conversion for types: (dtype(‘float64’), dtype(‘O’))
    cont_table = cont_table.astype('float64')
    genes = ca.row_coordinates(cont_table)

    gene_lengths = [len(record.seq) for record in filter_reference(parse(handle, 'fasta'), min_len_threshold, 'nuc')]

    x = genes.iloc[:, 0]
    y = genes.iloc[:, 1]
    fig = plt.figure(figsize=(16, 9))
    plt.scatter(x=x, y=y, alpha=0.8, c=gene_lengths, cmap='viridis', zorder=2)
    plt.grid(True, linestyle=':')
    plt.axvline(0, color='red', zorder=1)
    plt.axhline(0, color='red', zorder=1)
    plt.xlabel(f'Axis 0 (Variance: {round(ca.percentage_of_variance_[0], 4)}%)')
    plt.ylabel(f'Axis 1 (Variance: {round(ca.percentage_of_variance_[1], 4)}%)')
    c_bar = plt.colorbar()
    c_bar.set_label('Gene Length')
    plt.title(f'Total genes: {len(cont_table.iloc[:, 0])}')
    sup_title = f'CA of Codon RSCU of {organism_name} [genes]' if organism_name \
        else 'CA of Codon RSCU [genes]'
    plt.suptitle(sup_title)

    if save_image:
        make_dir(folder_path)
        name = f'CA_codon_rscu_gene_{organism_name}.png' if organism_name else 'CA_codon_rscu_gene.png'
        file_name = join(folder_path, name)
        if is_file_writeable(file_name):
            plt.savefig(file_name, dpi=500)

    plt.show()
    plt.close(fig)
