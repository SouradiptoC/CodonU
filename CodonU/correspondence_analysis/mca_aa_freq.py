import numpy as np
import pandas as pd
from prince import PCA
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.SeqIO import parse

from CodonU.analyzer.internal_comp import filter_reference


def mca_aa_freq(handle: str, genetic_table_num: int, min_len_threshold: int = 66, n_components: int = 20) -> \
        tuple[pd.DataFrame, np.ndarray]:
    """
    Calculates the contingency table and the inertia from amino acid frequency of gene

    **NOTE**: The fasta file must contain protein sequences in single letter format

    :param handle: Handle to the file, or the filename as a string
    :param genetic_table_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene (optional)
    :param n_components: The number of principal components to compute (optional)
    :return: The contingency table and inertia [inertia values lying between 0 and 1]
    """
    records = parse(handle, 'fasta')
    references = filter_reference(records, min_len_threshold)
    aas = [aa for aa, _ in unambiguous_dna_by_id[genetic_table_num].back_table.items()]
    aas.remove(None)  # as None is also returned to the back_table
    prot_names = [f'prot_{i}' for i in range(len(references))]
    contingency_table = pd.DataFrame(index=prot_names, columns=aas)
    for idx, prot in enumerate(references):
        for aa in aas:
            contingency_table[aa][prot_names[idx]] = prot.count(aa)
    pca = PCA(random_state=42, n_components=n_components)
    pca.fit(contingency_table)
    print('The inertia for respective components are:')
    for idx, inertia in enumerate(pca.explained_inertia_):
        print(f'Axis {idx + 1}: {inertia}')
    return contingency_table, pca.explained_inertia_
