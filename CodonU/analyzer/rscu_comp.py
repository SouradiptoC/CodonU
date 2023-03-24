from CAI import RSCU
from .internal_comp import filter_reference
from Bio.SeqIO import parse


def calculate_rscu(handle: str, genetic_code_num: int, min_len_threshold: int = 200, gene_analysis: bool = False) -> \
        dict[str, float | dict[str, float]]:
    """
    Calculates rscu values for each codon

    :param handle: Handle to the file, or the filename as a string
    :param genetic_code_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    :return: The dictionary containing codon and rscu value pairs if gene_analysis is false, otherwise the dictionary containing the gene name and the codon & rscu value pairs
    """
    records = parse(handle, 'fasta')
    references = filter_reference(records, min_len_threshold)
    if gene_analysis:
        rscu_dict = dict()
        for i, seq in enumerate(references):
            rscu_dict.update({f'gene_{i + 1}': RSCU([seq], genetic_code_num)})
        return rscu_dict
    else:
        reference = filter_reference(records, min_len_threshold)
        return RSCU(reference, genetic_code_num)
