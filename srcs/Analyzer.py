from statistics import median
from Bio.Seq import Seq
from Bio.SeqIO import parse
from Bio.Data.CodonTable import unambiguous_dna_by_id, NCBICodonTableDNA

from Bio.SeqUtils import GC123
from CAI import CAI, RSCU

from warnings import filterwarnings
from Errors import ThresholdError


def syn_codons(codon_table: NCBICodonTableDNA) -> dict[str, list[str]]:
    """
    Creates the protein, codon dictionary where protein is key

    :param codon_table: The codon table
    :return: The dict having protein as key
    """
    codon_dict = codon_table.forward_table
    syn_dict = dict()
    for codon, aa in codon_dict.items():
        syn_dict[aa] = syn_dict.get(aa, [])
        syn_dict[aa].append(codon)
    return syn_dict


def sf_vals(codon_table: NCBICodonTableDNA) -> dict[int, list[str]]:
    """
    Creates the sf value and protein dictionary where sf value is key

    :param codon_table: The codon table
    :return: The dict having sf values as key
    """
    syn_codon_dict = syn_codons(codon_table)
    sf_dic = dict()
    for aa, codons in syn_codon_dict.items():
        sf = len(codons)
        sf_dic[sf] = sf_dic.get(sf, [])
        sf_dic[sf].append(aa)
    return sf_dic


def cbi():
    pass


def gc_123(seq: Seq | str) -> tuple[float, float | int, float | int, float | int]:
    """
    Calculate G+C content: total, for first, second and third positions

    :param seq: Provided sequence
    :return: The tuple containing the result
    """
    return GC123(seq)


def filter_reference(records, threshold: float) -> list[Seq]:
    """
    Filters the list of reference based on given threshold of length

    :param records: A generator object holding the sequence objects
    :param threshold: Used to calculate minimum length required
    :return: The list of usable sequences
    """
    if not 0 <= threshold <= 1:
        raise ThresholdError
    reference = [record.seq for record in records]
    seq_len_lst = [len(seq) for seq in reference]
    min_len = median(seq_len_lst) * threshold
    filtered_lst = [seq for seq in reference if len(seq) >= min_len]
    return filtered_lst


def calculate_cai(records, genetic_code_num: int, threshold: float = 0.1) -> dict[str, float]:
    """
    Calculates cai for each codon

    :param records: The generator object containing sequence object
    :param genetic_code_num: Genetic table number for codon table
    :param threshold: Threshold value for filter
    :return: The dictionary containing codon and cai value pairs
    """
    reference = filter_reference(records, threshold)
    filterwarnings('ignore')
    cai_dict = dict()
    for codon in unambiguous_dna_by_id[genetic_code_num].forward_table:
        cai_val = CAI(codon, reference=reference, genetic_code=genetic_code_num)
        cai_dict.update({codon: cai_val})
    return cai_dict


def calculate_rscu(records, genetic_code_num: int, threshold: float = 0.1) -> dict[str, float]:
    """
    Calculates rscu values for each codon

    :param records: The generator object containing sequence object
    :param genetic_code_num: Genetic table number for codon table
    :param threshold: Threshold value for filter
    :return: The dictionary containing codon and rscu value pairs
    """
    reference = filter_reference(records, threshold)
    return RSCU(reference, genetic_code_num)


def calculate_cbi():
    pass


if __name__ == '__main__':
    print(gc_123('GGG'))
    handle = '/home/souro/Projects/final_yr/Results/Nucleotide/Staphylococcus_agnetis_nucleotide.fasta'
    lst = []
    records = parse(handle, 'fasta')
    # print(len(rscu(records, 1)))
    # for record in records:
    #     lst.append(record.seq)
    cai_dct = calculate_cai(records, 11)
    print(cai_dct['ATG'])
    # print(cai_dct.keys())
    # seq_lst = [record.seq for record in records]
    # filter_reference(records, 4.5)
