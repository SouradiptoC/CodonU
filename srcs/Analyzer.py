from statistics import median
from Bio.Seq import Seq
from Bio.SeqIO import parse
from Bio.Data.CodonTable import unambiguous_dna_by_id

from Bio.SeqUtils import GC123
from CAI import CAI, RSCU

from warnings import filterwarnings
from Errors import ThresholdError


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


def cai(records, genetic_code_num: int, threshold: float = 0.1) -> dict[str, float]:
    """
    Calculates cai for each codon

    :param records: The generator object containing sequence object
    :param genetic_code_num: Genetic table number for codon table
    :param threshold: Threshold value for filter
    :return:
    """
    reference = filter_reference(records, threshold)
    filterwarnings('ignore')
    cai_dict = dict()
    for codon in unambiguous_dna_by_id[genetic_code_num].forward_table:
        cai_val = CAI(codon, reference=reference, genetic_code=genetic_code_num)
        cai_dict.update({codon: cai_val})
    return cai_dict


def rscu(records, genetic_code_num: int, threshold: float = 0.1) -> dict[str, float]:
    """
    Calculates rscu values for each codon
    :param records: The generator object containing sequence object
    :param genetic_code_num: Genetic table number for codon table
    :param threshold: Threshold value for filter
    :return:
    """
    reference = filter_reference(records, threshold)
    return RSCU(reference, genetic_code_num)


#     6144 * 3456


if __name__ == '__main__':
    print(gc_123('GGG'))
    handle = '/home/souro/Projects/final_yr/Results/Nucleotide/Staphylococcus_agnetis_nucleotide.fasta'
    lst = []
    records = parse(handle, 'fasta')
    # print(len(rscu(records, 1)))
    # for record in records:
    #     lst.append(record.seq)
    cai_dct = cai(records, 11)
    print(cai_dct['ATG'])
    # print(cai_dct.keys())
    # seq_lst = [record.seq for record in records]
    # filter_reference(records, 4.5)
