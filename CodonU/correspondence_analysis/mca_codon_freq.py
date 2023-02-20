from collections import Counter
from itertools import chain
import pandas as pd
from prince import PCA
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.SeqIO import parse
from CodonU.analyzer.internal_comp import filter_reference
from CodonU.cua_warnings import NoCodonWarning


def mca_codon_freq(handle: str, genetic_code_num: int, min_len_threshold: int = 200, n_components: int = 59):
    records = parse(handle, 'fasta')
    references = filter_reference(records, min_len_threshold)
    codons = [codon for codon, _ in unambiguous_dna_by_id[genetic_code_num].forward_table.items()]
    gene_names = [f'gene_{i}' for i in range(len(references))]
    # len_lst = [len(gene) for gene in references]
    # max_len = max(len_lst)
    # s = list()
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
        # s.append(len(gene) / max_len * 100)
    pca = PCA(random_state=42, n_components=n_components)
    pca.fit(contingency_table)
    print('The inertias are:')
    for idx, inertia in enumerate(pca.explained_inertia_):
        print(f'Axis {idx + 1}: {inertia}')
