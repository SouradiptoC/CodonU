from .internal_comp import filter_reference, cbi
from warnings import filterwarnings
from Bio.Data.CodonTable import unambiguous_dna_by_id


def calculate_cbi(records, genetic_code_num: int, min_len_threshold: int = 200, gene_analysis: bool = False) -> \
        dict[str, tuple[float, str]] | dict[str, dict[str, tuple[float, str]]]:
    """
    Calculates cbi values for each amino acid

    :param records: The generator object containing sequence object
    :param genetic_code_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    :return: The dictionary containing amino acid and cbi value, optimal codon pairs
    """
    reference = filter_reference(records, min_len_threshold)
    filterwarnings('ignore')
    cbi_dict = dict()
    if gene_analysis:
        for i, seq in enumerate(reference):
            cbi_val_dict = dict()
            for codon in unambiguous_dna_by_id[genetic_code_num].forward_table:
                cbi_val = cbi(codon, reference=[seq], genetic_code=genetic_code_num)
                cbi_val_dict.update({codon: cbi_val})
            cbi_dict.update({f'gene_{i + 1}': cbi_val_dict})
    else:
        for aa in unambiguous_dna_by_id[genetic_code_num].protein_alphabet:
            cbi_val = cbi(aa, reference, genetic_code_num)
            cbi_dict.update({aa: cbi_val})
    return cbi_dict
