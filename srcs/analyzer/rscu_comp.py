from CAI import RSCU
from .internal_comp import filter_reference


def calculate_rscu(records, genetic_code_num: int, min_len_threshold: int = 200, gene_analysis: bool = False) -> \
        dict[str, float] | dict[str, dict[str, float]]:
    """
    Calculates rscu values for each codon

    :param records: The generator object containing sequence object
    :param genetic_code_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    :return: The dictionary containing codon and rscu value pairs
    """
    references = filter_reference(records, min_len_threshold)
    if gene_analysis:
        rscu_dict = dict()
        for i, seq in enumerate(references):
            rscu_dict.update({f'gene_{i + 1}': RSCU([seq], genetic_code_num)})
        return rscu_dict
    else:
        reference = filter_reference(records, min_len_threshold)
        return RSCU(reference, genetic_code_num)
