import pandas as pd
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.SeqIO import parse
from CodonU.analyzer import filter_reference
from CodonU.file_handler import make_dir
from CodonU.file_handler.internal_comp import is_file_writeable
from collections import Counter
from itertools import chain
from os.path import join, abspath


def build_contingency_table_codon_count(handle: str, genetic_table_num: int, min_len_threshold: int = 200,
                                        save_file: bool = False, file_name: str = 'contingency_codon_count',
                                        folder_path: str = 'Report') -> pd.DataFrame:
    """
    Creates the contingency table of codon frequency\n
    **Note** Gene descriptions are index headings and codons are column headings

    :param handle: Handle to the file, or the filename as a string
    :param genetic_table_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :param save_file: Option for saving the values in xlsx format (Optional)
    :param file_name: Intended file name (Optional)
    :param folder_path: Folder path where image should be saved (optional)
    :return: The contingency table as a pandas DataFrame object
    """
    records = parse(handle, 'fasta')
    filtered_records = filter_reference(records, min_len_threshold=min_len_threshold, _type='nuc')
    codons = list(unambiguous_dna_by_id[genetic_table_num].forward_table.keys())
    genes = [record.description for record in filtered_records]

    contingency_table = pd.DataFrame(index=genes, columns=codons)
    for idx, gene in enumerate(filtered_records):
        sequences = ((seq[i:i + 3].upper() for i in range(0, len(seq), 3)) for seq in [gene.seq])
        _codons = chain.from_iterable(sequences)
        counts = Counter(_codons)
        for codon in codons:
            contingency_table[codon][genes[idx]] = counts[codon]

    if save_file:
        name = file_name + '.xlsx'
        make_dir(folder_path)
        file_path = join(folder_path, name)
        if is_file_writeable(file_path):
            contingency_table.to_excel(file_path, float_format='%.5f', columns=contingency_table.columns)
        print(f'The Contingency table of codon count can be found at: {abspath(file_path)}')

    return contingency_table
