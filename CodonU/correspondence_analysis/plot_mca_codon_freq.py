from collections import Counter
from itertools import chain
from os.path import join
import pandas as pd
import matplotlib.pyplot as plt
from prince import PCA
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.SeqIO import parse
from CodonU.analyzer.internal_comp import filter_reference
from CodonU.cua_warnings import NoCodonWarning


def plot_mca_codon_freq(handle: str, genetic_table_num: int, min_len_threshold: int = 200, n_components: int = 59,
                        organism_name: str | None = None, save_image: bool = False, folder_path: str = ''):
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
    codons = [codon for codon, _ in unambiguous_dna_by_id[genetic_table_num].forward_table.items()]
    gene_names = [f'gene_{i}' for i in range(len(references))]
    len_lst = [len(gene) for gene in references]
    max_len = max(len_lst)
    s = list()
    contingency_table = pd.DataFrame(index=gene_names, columns=codons)
    for idx, gene in enumerate(references):
        sequences = ((sequence[i:i + 3].upper() for i in range(0, len(sequence), 3)) for sequence in [gene])
        _codons = chain.from_iterable(sequences)
        counts = Counter(_codons)
        for codon in codons:
            try:
                contingency_table[codon][gene_names[idx]] = counts[codon]
            except KeyError:
                warn = NoCodonWarning(codon)
                warn.warn()
                contingency_table[codon][gene_names[idx]] = 0
        s.append(len(gene) / max_len * 100)
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
        name = f'Multivariate_analysis_codon_freq_{organism_name}.png' if organism_name else 'Multivariate_analysis_codon_freq.png'
        file_name = join(folder_path, name)
        plt.savefig(file_name, dpi=500)
    plt.show()
    plt.close()
