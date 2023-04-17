from cai2 import CAI
from warnings import filterwarnings
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.SeqIO import parse
from .internal_comp import filter_reference


def calculate_cai(handle: str, genetic_code_num: int, min_len_threshold: int = 200, gene_analysis: bool = False) -> \
        dict[str, float | dict[str, float]]:
    """
    Calculates cai values for each codon

    :param handle: Handle to the file, or the filename as a string
    :param genetic_code_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    :return: The dictionary containing codon and cai value pairs if gene_analysis is False, otherwise returns the
    dictionary containing gene name and corresponding codon and cai value pairs
    """
    filterwarnings('ignore')
    cai_dict = dict()
    records = parse(handle, 'fasta')
    reference = filter_reference(records, min_len_threshold)
    if gene_analysis:
        for i, seq in enumerate(reference):
            cai_val_dict = dict()
            for codon in unambiguous_dna_by_id[genetic_code_num].forward_table:
                cai_val = CAI(codon, reference=[seq], genetic_code=genetic_code_num)
                cai_val_dict.update({codon: cai_val})
            cai_dict.update({f'gene_{i + 1}': cai_val_dict})
    else:
        for codon in unambiguous_dna_by_id[genetic_code_num].forward_table:
            cai_val = CAI(codon, reference=reference, genetic_code=genetic_code_num)
            cai_dict.update({codon: cai_val})
    return cai_dict
