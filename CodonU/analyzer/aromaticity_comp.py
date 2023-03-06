from warnings import filterwarnings
from CodonU.analyzer.internal_comp import filter_reference, aromaticity


def calculate_aromaticity(records, min_len_threshold: int = 66, gene_analysis: bool = False) -> \
        dict[str, float] | float:
    """
    Calculates the aromaticity score for a given protein sequence

    :param records: The generator object containing protein sequence object
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    :return: The aromaticity score of given sequence if gene_analysis is false, else the dictionary containing gene number and corresponding GRAVY score
    """
    references = filter_reference(records, min_len_threshold)
    filterwarnings('ignore')
    if gene_analysis:
        gravy_dict = dict()
        for i, seq in enumerate(references):
            gravy_dict.update({f'prot_seq{i + 1}': aromaticity([seq])})
        return gravy_dict
    else:
        return aromaticity(references)
