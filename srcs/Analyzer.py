from math import nan
from statistics import median
from Bio.Seq import Seq
from Bio.SeqIO import parse
from Bio.Data.CodonTable import unambiguous_dna_by_id, NCBICodonTableDNA

from Bio.SeqUtils import GC123, seq3
from CAI import CAI, RSCU
from itertools import chain
from collections import Counter

from warnings import filterwarnings
from Errors import ThresholdError, MissingCodonWarning, NoSynonymousCodonWarning


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


def cbi(prot_seq: Seq | str, reference: list[Seq], genetic_code: int) -> tuple[float, str]:
    """
    Calculates codon bias index (CBI) for a given protein seq

    :param prot_seq: The Amino Acid
    :param reference: List of reference nucleotide sequences
    :param genetic_code: Genetic table number for codon table
    :return: A tuple of CBI val and the optimal codon
    :raises NoSynonymousCodonWarning: When there is no synonymous codons
    :raises MissingCodonWarning: When no codons translate to provided Amino acid
    """
    sequences = ((sequence[i:i + 3].upper() for i in range(0, len(sequence), 3)) for sequence in reference)
    codons = chain.from_iterable(sequences)
    counts = Counter(codons)
    syn_codon_dict = syn_codons(unambiguous_dna_by_id[genetic_code])
    sf_val_dict = sf_vals(unambiguous_dna_by_id[genetic_code])
    sorted_keys = sorted(sf_val_dict.keys(), reverse=True)
    sorted_sf_dict = {key: sf_val_dict[key] for key in sorted_keys}
    cbi_val, opt_codon = nan, None
    for num, aa_lst in sorted_sf_dict.items():
        if num == 1:
            opt_codon = syn_codon_dict[prot_seq][0]
            warn = NoSynonymousCodonWarning(seq3(prot_seq))
            warn.warn()
            break
        if prot_seq in aa_lst:
            codon_lst = syn_codon_dict[prot_seq]
            count_lst = [counts[codon] for codon in codon_lst]
            tot_count = sum(count_lst)
            ran_count = tot_count / num
            opt_count = max(count_lst)
            try:
                cbi_val = (opt_count - ran_count) / (tot_count - ran_count)
                opt_codon = codon_lst[count_lst.index(max(count_lst))]
            except ZeroDivisionError:
                warn = MissingCodonWarning(seq3(prot_seq))
                warn.warn()
            break

    return cbi_val, opt_codon


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


def calculate_cbi(records, genetic_code_num: int, threshold: float = 0.1) -> dict[str, tuple[float, str]]:
    """
    Calculates cbi values for each amino acid

    :param records: The generator object containing sequence object
    :param genetic_code_num: Genetic table number for codon table
    :param threshold: Threshold value for filter
    :return: The dictionary containing amino acid and cbi value, optimal codon pairs
    """
    reference = filter_reference(records, threshold)
    # filterwarnings('ignore')
    cbi_dict = dict()
    for aa in unambiguous_dna_by_id[genetic_code_num].protein_alphabet:
        cbi_val = cbi(aa, reference, genetic_code_num)
        cbi_dict.update({aa: cbi_val})
    return cbi_dict


if __name__ == '__main__':
    # handle = '/home/souro/Projects/final_yr/Results/Nucleotide/Staphylococcus_agnetis_nucleotide.fasta'
    handle = '/home/souro/Projects/final_yr/Results/Nucleotide/temp.fasta'
    # records = parse(handle, 'fasta')
    # cai_dict = calculate_cai(records, 11)
    # for i, j in cai_dict.items():
    #     print(i, j)
    print('#####################################################')
    records = parse(handle, 'fasta')
    cbi_dict = calculate_cbi(records, 11)
    for i, j in cbi_dict.items():
        print(seq3(i), j)
