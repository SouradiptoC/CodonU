from os.path import join
import matplotlib.pyplot as plt
from prince import PCA
from Bio.SeqIO import parse
from CodonU.analyzer.internal_comp import filter_reference
from CodonU.correspondence_analysis.mca_codon_freq import mca_codon_freq
from CodonU.file_handler import make_dir
from CodonU.file_handler.internal_comp import is_file_writeable


def plot_mca_codon_freq(handle: str, genetic_table_num: int, min_len_threshold: int = 200, n_components: int = 59,
                        organism_name: str | None = None, save_image: bool = False, folder_path: str = 'Report'):
    """
    Plots the principal component analysis based on codon frequency

    :param handle: Handle to the file, or the filename as a string
    :param genetic_table_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene (optional)
    :param n_components: The number of principal components to compute (optional)
    :param organism_name: Name of organism (optional)
    :param save_image: Options for saving the image (optional)
    :param folder_path: Folder path where image should be saved (optional)
    """
    records = parse(handle, 'fasta')
    references = filter_reference(records, min_len_threshold)
    len_lst = [len(gene) for gene in references]
    max_len = max(len_lst)
    s = [len(gene) / max_len * 100 for gene in references]
    contingency_table, _ = mca_codon_freq(handle, genetic_table_num, min_len_threshold, n_components)
    pca = PCA(random_state=42, n_components=n_components)
    pca.fit(contingency_table)
    plot_df = pca.row_coordinates(contingency_table)
    x = plot_df.iloc[:, 0]
    y = plot_df.iloc[:, 1]
    plt.figure(figsize=(9, 5.25))
    plt.scatter(x, y, s, alpha=0.5, c=len_lst, cmap='viridis', zorder=2)
    plt.grid(True, linestyle=':')
    plt.axvline(0, color='red', zorder=1)
    plt.axhline(0, color='red', zorder=1)
    plt.xlabel(f'Axis 0 (inertia: {round(pca.explained_inertia_[0] * 100, 4)}%)')
    plt.ylabel(f'Axis 1 (inertia: {round(pca.explained_inertia_[1] * 100, 4)}%)')
    c_bar = plt.colorbar()
    c_bar.set_label('Length of gene')
    plt.title(f'Total genes: {len(references)}')
    sup_title = f'Multivariate analysis of Codon Frequency of {organism_name}' if organism_name else 'Multivariate analysis of Codon Frequency'
    plt.suptitle(sup_title)
    if save_image:
        make_dir(folder_path)
        name = f'Multivariate_analysis_codon_freq_{organism_name}.png' if organism_name else 'Multivariate_analysis_codon_freq.png'
        file_name = join(folder_path, name)
        if is_file_writeable(file_name):
            plt.savefig(file_name, dpi=500)
    plt.show()
    plt.close()
