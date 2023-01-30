from math import isnan
from statistics import median
from Bio.Seq import Seq
from Bio.SeqIO import parse
from Bio.Data.IUPACData import unambiguous_dna_letters

from Bio.SeqUtils import GC, GC123
from CAI import CAI

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
    reference = list()
    for record in records:
        reference.append(record.seq)
    seq_len_lst = [len(seq) for seq in reference]
    min_len = median(seq_len_lst) * threshold
    filtered_lst = [seq for seq in reference if len(seq) >= min_len]
    return filtered_lst


def cai(records, genetic_code: int, threshold: float = 0.1) -> dict[str, float]:
    """

    :param records:
    :param genetic_code:
    :param threshold:
    :return:
    """
    reference = filter_reference(records, threshold)
    filterwarnings('ignore')
    cai_dict = dict()
    for b1 in unambiguous_dna_letters:
        for b2 in unambiguous_dna_letters:
            for b3 in unambiguous_dna_letters:
                codon = b1 + b2 + b3
                cai_val = CAI(codon, reference=reference, genetic_code=genetic_code)
                if not isnan(cai_val):
                    cai_dict.update({codon: cai_val})
    return cai_dict


#     6144 * 3456


if __name__ == '__main__':
    print(gc_123('GGG'))
    handle = '/home/souro/Projects/final_yr/Results/Nucleotide/Staphylococcus_agnetis_nucleotide.fasta'
    lst = []
    # records = parse(handle, 'fasta')
    # for record in records:
    #     lst.append(record.seq)
    # cai_dct = cai(records, 11)
    # print(cai_dct.keys())
    filter_reference(records, 4.5)
