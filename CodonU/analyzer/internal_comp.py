from collections import Counter
from itertools import chain
from math import nan, isnan
from statistics import mean
from Bio.Data.CodonTable import NCBICodonTableDNA, unambiguous_dna_by_id
from Bio.Seq import Seq
from Bio.SeqUtils import seq3, GC123
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from CodonU.cua_warnings import NoSynonymousCodonWarning, MissingCodonWarning
from CodonU.cua_errors import NoProteinError


def g3(seq: Seq | str) -> float:
    """
    Calculates percentage of G content for third position

    :param seq: Provided sequence
    :return: Percentage of G content
    """
    g3_val = 0
    for i in range(0, len(seq), 3):
        codon = seq[i: i + 3]
        if len(codon) < 3:
            codon += '  '
        if codon[-1] in 'Gg':
            g3_val += 1
    return g3_val / (len(seq) / 3) * 100


def a3(seq: Seq | str) -> float:
    """
    Calculates percentage of A content for third position

    :param seq: Provided sequence
    :return: Percentage of A content
    """
    a3_val = 0
    for i in range(0, len(seq), 3):
        codon = seq[i: i + 3]
        if len(codon) < 3:
            codon += '  '
        if codon[-1] in 'Aa':
            a3_val += 1
    return a3_val / (len(seq) / 3) * 100


def gc_123(seq: Seq | str) -> tuple[float, float | int, float | int, float | int]:
    """
    Calculate G+C content: total, for first, second and third positions

    :param seq: Provided sequence
    :return: The G+C percentage for the entire sequence, and the three codon positions
    """
    return GC123(seq)


def at_123(seq: Seq | str) -> tuple[float, float | int, float | int, float | int]:
    """
    Calculate G+C content: total, for first, second and third positions

    :param seq: Provided sequence
    :return: The A+T percentage for the entire sequence, and the three codon positions
    """
    gc_tot, gc_1, gc_2, gc_3 = gc_123(seq)
    return 100 - gc_tot, 100 - gc_1, 100 - gc_2, 100 - gc_3


def filter_reference(records, min_len_threshold: int) -> list[Seq]:
    """
    Filters the list of reference based on given threshold of length

    :param records: A generator object holding the sequence objects
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :return: The list of usable sequences
    """
    reference = [record.seq for record in records]
    filtered_lst = [seq for seq in reference if len(seq) >= min_len_threshold]
    return filtered_lst


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
    sorted_keys = sorted(sf_dic.keys(), reverse=True)
    sorted_sf_dict = {key: sf_dic[key] for key in sorted_keys}
    return sorted_sf_dict


def cbi(prot_seq: Seq | str, reference: list[Seq], genetic_code: int) -> tuple[float, str]:
    """
    Calculates codon bias index (CBI) for a given protein seq based on Bennetzen and Hall (1982)

    :param prot_seq: The Protein Sequence
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
    cbi_val, opt_codon = nan, None
    for num, aa_lst in sf_val_dict.items():
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


def enc(references: list[Seq | str], genetic_code: int) -> float:
    """
    Calculates Effective number of codons (Enc) based on Wright (1989) and Fuglsang (2004)

    :param references: List of reference nucleotide sequences
    :param genetic_code: Genetic table number for codon table
    :return: Calculated Enc value for the sequence(s)
    :raises MissingCodonWarning: If there is no codon for a certain amino acid
    :raises NoProteinError: If there is no codon for a certain set of amino acid
    """
    sequences = ((sequence[i: i + 3].upper() for i in range(0, len(sequence), 3)) for sequence in references)
    codons = chain.from_iterable(sequences)
    counts = Counter(codons)
    syn_dct = syn_codons(unambiguous_dna_by_id[genetic_code])
    sf_dct = sf_vals(unambiguous_dna_by_id[genetic_code])
    F_val_dict = dict()
    for num, aa_lst in sf_dct.items():
        F_val_dict.update({num: [0.0 for _ in range(len(aa_lst))]})
    # [sf_6 avg, sf_4 avg, sf_3 avg, sf_2 avg, sf_1 avg] for F_val_avg_lst
    F_val_avg_lst = list()
    for num, aa_lst in sf_dct.items():
        if num == 1:
            F_val_dict.update({num: [1.0, 1.0]})
            break
        for aa in aa_lst:
            codon_lst = syn_dct[aa]
            count_lst = [counts[codon] for codon in codon_lst]
            tot_count = float(sum(count_lst))
            p_2 = 0.0
            F = 0.0
            for count in count_lst:
                try:
                    p_2 += (count / tot_count) ** 2
                except ZeroDivisionError:
                    warn = MissingCodonWarning(seq3(aa))
                    warn.warn()
                    F = nan
                    break
            if not isnan(F):
                try:
                    F = ((tot_count * p_2) - 1) / (tot_count - 1)
                except ZeroDivisionError:
                    F = 1.0  # lim(x->0) x/x = 1 and tot_count = 0
            F_val_lst = F_val_dict.get(num)
            F_val_lst[aa_lst.index(aa)] = F
            F_val_dict.update({num: F_val_lst})
    for num, val_lst in F_val_dict.items():
        if nan not in val_lst and num != 3:
            F_val_avg_lst.append(mean(val_lst))
        elif nan in val_lst and num != 3:
            refined_lst = [val for val in val_lst if val is not nan]
            F_val_avg_lst.append(mean(refined_lst))
        else:
            F_val_avg_lst.append(nan)
    if F_val_avg_lst[2] in [0, nan]:
        f_2, f_4, f_6 = F_val_avg_lst[3], F_val_avg_lst[1], F_val_avg_lst[0]
        if 0 not in [f_2, f_4, f_6]:
            f_3 = ((((2 / f_2) - 1) ** -1) + (((2 / (3 * f_4)) + (1 / 3)) ** -1) + (
                    ((2 / (5 * f_6)) + (3 / 5)) ** -1)) / 3.0
        else:
            raise NoProteinError(references[0])
        F_val_avg_lst[2] = f_3
    # [sf_6 avg, sf_4 avg, sf_3 avg, sf_2 avg, sf_1 avg]
    enc_val = 2 + (9 / F_val_avg_lst[3]) + (1 / F_val_avg_lst[2]) + (5 / F_val_avg_lst[1]) + (3 / F_val_avg_lst[0])
    return enc_val if enc_val < 61 else 61.00


def gravy(seq: Seq | str) -> float:
    """
    Computes the GRAVY score according to Kyte and Doolittle (1982)

    :param seq: Protein sequence
    :return: The GRAVY score
    """
    if not isinstance(seq, str):
        _seq = str(seq)
        seq = _seq
    return ProteinAnalysis(seq).gravy()


def aromaticity(seq: Seq | str) -> float:
    """
    Calculate the aromaticity score according to Lobry (1994).

    :param seq: Protein sequence
    :return: The aromaticity score
    """
    if not isinstance(seq, str):
        _seq = str(seq)
        seq = _seq
    return ProteinAnalysis(seq).aromaticity()
