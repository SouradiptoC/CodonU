from collections import Counter
from itertools import chain
import pandas as pd
import numpy as np
from prince import PCA
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.SeqIO import parse
from CodonU.analyzer.internal_comp import filter_reference
from CodonU.cua_warnings import NoCodonWarning


def mca_codon_freq(handle: str, genetic_table_num: int, min_len_threshold: int = 200, n_components: int = 59) -> \
        tuple[pd.DataFrame, np.ndarray]:
    """
    Calculates the contingency table and the inertia from codon frequency of gene

    :param handle: Handle to the file, or the filename as a string
    :param genetic_table_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene (optional)
    :param n_components: The number of principal components to compute (optional)
    :return: The contingency table and inertia [inertia values lying between 0 and 1]
    """
    records = parse(handle, 'fasta')
    references = filter_reference(records, min_len_threshold)
    codons = [codon for codon, _ in unambiguous_dna_by_id[genetic_table_num].forward_table.items()]
    gene_names = [f'gene_{i}' for i in range(len(references))]
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
    pca = PCA(random_state=42, n_components=n_components)
    pca.fit(contingency_table)
    print('The inertia for respective components are:')
    for idx, inertia in enumerate(pca.explained_inertia_):
        print(f'Axis {idx + 1}: {inertia}')
    return contingency_table, pca.explained_inertia_
