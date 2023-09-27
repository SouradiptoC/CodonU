from collections import Counter
from itertools import chain
from math import nan, isnan
from statistics import mean
from typing import Optional
from Bio.Data.CodonTable import NCBICodonTableDNA, unambiguous_dna_by_id, unambiguous_dna_by_name, register_ncbi_table
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq3, GC123
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.stats import gmean
from CodonU.cua_warnings import NoSynonymousCodonWarning, MissingCodonWarning
from CodonU.cua_errors import NoProteinError, CodonTableExistsError, BadSequenceError, NucleotideError


def is_not_bad_seq(seq: Seq | str, code: int, _type: str) -> bool:
    """
    Checks if the sequence is bad i.e. length of the sequence is not divisible by 3

    :param seq: The nucleotide sequence
    :param code: The code to call BadSequenceError (1 or 2)
    :param _type: Type of sequence, i.e. 'nuc'
    :return: True if seq is not bad
    :raises BadSequenceError: If the seq is bad
    """
    if _type == 'nuc':
        if len(seq) % 3 == 0:
            return True
        else:
            raise BadSequenceError(seq, code)
    elif _type == 'aa':
        return True


def not_contains_amb_letter(seq: Seq | str) -> bool:
    """
    Checks if provided sequence contains ambiguous DNA letters

    :param seq: Provided sequence
    :return: True if sequence does not contain ambiguous letter
    :raise NucleotideError: If sequence contain ambiguous letter
    """
    amb_letters = "BDHKMNRVXY"
    if any(letter in amb_letters for letter in seq):
        raise NucleotideError(1)
    return True


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


def custom_codon_table(name: str, alt_name: Optional[str], genetic_code_id: int, forward_table: dict[str, str],
                       start_codons: list[str], stop_codons: list[str]) -> None:
    """
    Registers a new Codon Table as provided by the user. \n
    **Note**: The scope of the newly registered table is limited to the working file only

    :param name: Name for the table
    :param alt_name: Short name for the table
    :param genetic_code_id: Genetic code number for the table
    :param forward_table: A dict containing mapping of codons to proteins [excluding stop codons]
    :param start_codons: A list of possible start codons
    :param stop_codons: A list of possible stop codons
    :raises CodonTableExistsError: If the name, alt_name or genetic_code_id already exists
    """
    if genetic_code_id in unambiguous_dna_by_id.keys():
        raise CodonTableExistsError(1, genetic_code_id)
    if name in unambiguous_dna_by_name.keys():
        raise CodonTableExistsError(2, name)
    if alt_name in unambiguous_dna_by_name.keys():
        raise CodonTableExistsError(3, alt_name)

    register_ncbi_table(name=name, alt_name=alt_name, id=genetic_code_id, table=forward_table,
                        start_codons=start_codons, stop_codons=stop_codons)


def filter_reference(records, min_len_threshold: int, _type: str) -> list[SeqRecord]:
    """
    Filters the list of reference based on given threshold of length

    :param records: A generator object holding the sequence objects
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :param _type: Type of sequence, i.e. 'nuc' or 'aa
    :return: The list of usable sequences
    """
    filtered_lst = [record for record in records if
                    len(record.seq) >= min_len_threshold and is_not_bad_seq(record, 2, _type)]
    return filtered_lst


def reverse_table(codon_table: NCBICodonTableDNA) -> dict[str, list[str]]:
    """
    Creates the protein, codon dictionary where protein is key \n
    e.g. 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']

    :param codon_table: The codon table
    :return: The dict having protein as key
    """
    codon_dict = codon_table.forward_table
    reverse_table_dict = dict()
    for codon, aa in codon_dict.items():
        reverse_table_dict[aa] = reverse_table_dict.get(aa, [])
        reverse_table_dict[aa].append(codon)
    return reverse_table_dict


def syn_codons(codon_table: NCBICodonTableDNA) -> dict[str, list[str]]:
    """
    Creates the codon, synonymous codon family dictionary where codon is the key \n
    e.g. 'TTA': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']

    :param codon_table: The codon table
    :return: The dict having individual codons as keys
    """
    rev_table = reverse_table(codon_table)
    return {codon: rev_table[codon_table.forward_table[codon]] for codon in codon_table.forward_table.keys()}


def sf_vals(codon_table: NCBICodonTableDNA) -> dict[int, list[str]]:
    """
    Creates the sf value and protein dictionary where sf value is key \n
    e.g. 6: ['L', 'S', 'R']

    :param codon_table: The codon table
    :return: The dict having sf values as key
    """
    syn_codon_dict = reverse_table(codon_table)
    sf_dic = dict()
    for aa, codons in syn_codon_dict.items():
        sf = len(codons)
        sf_dic[sf] = sf_dic.get(sf, [])
        sf_dic[sf].append(aa)
    sorted_keys = sorted(sf_dic.keys(), reverse=True)
    sorted_sf_dict = {key: sf_dic[key] for key in sorted_keys}
    return sorted_sf_dict


def rscu(references: list[Seq | str], genetic_code: int) -> dict[str, float]:
    """
    Calculates relative synonymous codon usage (RSCU) value for a given nucleotide sequence according to Sharp and Li (1987)

    :param references: List of reference nucleotide sequences
    :param genetic_code: Genetic table number for codon table
    :return: A dictionary containing codons and their respective RSCU values
    """
    sequences = ((sequence[i:i + 3].upper() for i in range(0, len(sequence), 3)) for sequence in references if
                 is_not_bad_seq(sequence, 1, 'nuc'))
    codons = chain.from_iterable(sequences)  # flat list for Counter
    counts = dict(Counter(codons))  # Counters can not handle float values, hence dict

    # "Note that if a certain codon is never used in the reference set then the CAI for any other gene in which that
    # codon appears becomes zero. To overcome this problem
    # we assign a value of 0.5 to any X_{ij} that would otherwise be zero." Sharp and Li (1987) [pp. 1285]

    for codon in unambiguous_dna_by_id[genetic_code].forward_table:
        if codon not in counts.keys():
            counts[codon] = 0.5

    _syn_codons = syn_codons(unambiguous_dna_by_id[genetic_code])

    rscu_dict = dict()
    for codon in unambiguous_dna_by_id[genetic_code].forward_table:
        rscu_dict[codon] = counts[codon] / (
                (len(_syn_codons[codon]) ** -1) * (sum((counts[_codon] for _codon in _syn_codons[codon])))
        )

    return rscu_dict


def weights_for_cai(references: list[Seq | str], genetic_code: int) -> dict[str, float]:
    """
    Calculates relative adaptiveness/weight value for a given nucleotide sequence according to Sharp and Li (1987)

    :param references: List of reference nucleotide sequences
    :param genetic_code: Genetic table number for codon table
    :return: A dictionary containing codons and their respective weights
    """
    rscu_dict = rscu(references, genetic_code)
    _syn_codons = syn_codons(unambiguous_dna_by_id[genetic_code])

    weights = dict()
    for codon in rscu_dict:
        weights[codon] = rscu_dict[codon] / max(rscu_dict[_codon] for _codon in _syn_codons[codon])

    return weights


def cai(nuc_seq: Seq | str, references: list[Seq | str], genetic_code: int) -> float:
    """
    Calculates Codon Adaptive Index (CAI) value for a given nucleotide sequence according to Sharp and Li (1987)

    :param nuc_seq: The Nucleotide Sequence
    :param references: List of reference nucleotide sequences
    :param genetic_code: Genetic table number for codon table
    :return: The CAI value for given sequence
    """
    sequences = [nuc_seq[i: i + 3].upper() for i in range(0, len(nuc_seq), 3) if is_not_bad_seq(nuc_seq, 1, 'nuc')]
    weights_dict = weights_for_cai(references, genetic_code)
    seq_weights = list()
    non_syn_codons = [codon for codon in syn_codons(unambiguous_dna_by_id[genetic_code]).keys() if
                      len(syn_codons(unambiguous_dna_by_id[genetic_code])[codon]) == 1]

    for codon in sequences:
        if codon not in non_syn_codons:
            seq_weights.append(weights_dict[codon])

    return float(gmean(seq_weights))


def cbi(prot_seq: Seq | str, references: list[Seq | str], genetic_code: int) -> tuple[float, str]:
    """
    Calculates codon bias index (CBI) for a given protein seq based on Bennetzen and Hall (1982)

    :param prot_seq: The Protein Sequence
    :param references: List of reference nucleotide sequences
    :param genetic_code: Genetic table number for codon table
    :return: A tuple of CBI val and the optimal codon
    :raises NoSynonymousCodonWarning: When there is no synonymous codons
    :raises MissingCodonWarning: When no codons translate to provided Amino acid
    """
    sequences = ((sequence[i:i + 3].upper() for i in range(0, len(sequence), 3)) for sequence in references if
                 is_not_bad_seq(sequence, 1, 'nuc'))
    codons = chain.from_iterable(sequences)
    counts = Counter(codons)
    syn_codon_dict = reverse_table(unambiguous_dna_by_id[genetic_code])
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
    sequences = ((sequence[i: i + 3].upper() for i in range(0, len(sequence), 3)) for sequence in references if
                 is_not_bad_seq(sequence, 1, 'nuc'))
    codons = chain.from_iterable(sequences)
    counts = Counter(codons)
    syn_dct = reverse_table(unambiguous_dna_by_id[genetic_code])
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
        # correction based on Fuglsang (2004)
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
