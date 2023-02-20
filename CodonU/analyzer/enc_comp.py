from warnings import filterwarnings
from .internal_comp import filter_reference, enc


def calculate_enc(records, genetic_code_num: int, min_len_threshold=200, gene_analysis: bool = False) -> \
        float | dict[str, float]:
    """
    Calculates ENc value for a given sequences

    :param records: The generator object containing sequence object
    :param genetic_code_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    :return: The ENc value or a dictionary containing gene number and corresponding ENc value
    """
    references = filter_reference(records, min_len_threshold)
    filterwarnings('ignore')
    if gene_analysis:
        enc_dict = dict()
        for i, seq in enumerate(references):
            enc_dict.update({f'gene_{i + 1}': enc([seq], genetic_code_num)})
        return enc_dict
    else:
        return enc(references, genetic_code_num)
