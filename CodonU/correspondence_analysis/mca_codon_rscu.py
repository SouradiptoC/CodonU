import numpy as np
import pandas as pd
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.SeqIO import parse
from prince import PCA
from CodonU.analyzer import calculate_rscu


def mca_codon_rscu(handle: str, genetic_table_num: int, min_len_threshold: int = 200, n_components: int = 59) -> \
        tuple[pd.DataFrame, np.ndarray]:
    """
    Calculates the contingency table and the inertia from RSCU of every codon of every genes

    :param handle: Handle to the file, or the filename as a string
    :param genetic_table_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene (optional)
    :param n_components: The number of principal components to compute (optional)
    :return: The contingency table and inertia [inertia values lying between 0 and 1]
    """
    records = parse(handle, 'fasta')
    codons = [codon for codon, _ in unambiguous_dna_by_id[genetic_table_num].forward_table.items()]
    rscu_dict = calculate_rscu(records, genetic_table_num, min_len_threshold, gene_analysis=True)
    gene_names = list(rscu_dict.keys())
    contingency_table = pd.DataFrame(index=gene_names, columns=codons)
    for gene in gene_names:
        for codon in codons:
            contingency_table[codon][gene] = rscu_dict[gene][codon]
    pca = PCA(random_state=42, n_components=n_components)
    pca.fit(contingency_table)
    print('The inertia for respective components are:')
    for idx, inertia in enumerate(pca.explained_inertia_):
        print(f'Axis {idx + 1}: {inertia}')
    return contingency_table, pca.explained_inertia_
