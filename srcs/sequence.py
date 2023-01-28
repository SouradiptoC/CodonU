from Bio.Seq import Seq
import Bio.Data.CodonTable as cd_tab
import Bio.Seq
from Bio.Data.IUPACData import unambiguous_dna_letters, unambiguous_rna_letters, ambiguous_dna_letters, \
    ambiguous_rna_letters
from Errors import InternalStopCodonError, NucleotideError, TerSeqError
import pandas as pd


def validate_internal_stop_codon(genome: str | Bio.Seq.Seq, table_num: int):
    """
    Validates if internal stop codons are present
    :param genome: The sequence
    :param table_num: Table number
    :return:
    """
    for i in range(0, len(genome) - 3, 3):
        codon = genome[i: i + 3]
        if codon in cd_tab.generic_by_id[table_num].stop_codons:
            raise InternalStopCodonError(i)
    return True


def validate_nucleotide(genome: str | Bio.Seq.Seq, rna_flag=False):
    """
    Validates Nucleotides
    :param genome: The sequence
    :param rna_flag: flag for RNA
    :return:
    """
    bases = unambiguous_rna_letters if rna_flag else unambiguous_dna_letters
    amb_bases = ambiguous_rna_letters if rna_flag else ambiguous_dna_letters
    for i in range(len(genome)):
        if genome[i] not in bases:
            if genome[i] in amb_bases:
                raise NucleotideError(1, i)
            else:
                raise NucleotideError(2, i)
    return True


def validate_terminal_codon(genome: str | Bio.Seq.Seq, table_num: int):
    """
    Validates terminal codons viz start and stop codons
    :param genome: The sequence
    :param table_num: Table number
    :return:
    """
    if genome[:3] not in cd_tab.generic_by_id[table_num].start_codons:
        raise TerSeqError(1)
    elif genome[-3:] not in cd_tab.generic_by_id[table_num].stop_codons:
        raise TerSeqError(-1)
    return True


class Sequence:
    def __init__(self, seq, table_num):
        # TODO Check if already seq type or not
        # if validate_nucleotide(seq) and validate_terminal_codon(seq, table_num) and validate_internal_stop_codon(seq,
        #                                                                                                          table_num):
        self.seq = Seq(seq)

    def tot_kmer(self, k_len: int = 3):
        return len(self.seq) // k_len

    def count_gc(self):
        return (self.seq.count('G') + self.seq.count('C')) / len(self.seq)

    def count_at(self):
        return 1 - self.count_gc()

    def count_gc3(self):
        gc_str = 'GC'
        count, k, k_len = 0, 0, 3
        for _ in range(self.tot_kmer()):
            _seq = self.seq[k: k + k_len]
            k += k_len
            if _seq[2] in gc_str:
                count += 1
        return count / self.tot_kmer()

    def count_gc_not_gc3(self):
        pass

    # TODO implement count for gc_not_gc3 and individual counts for 3 positions

    def position_prob(self):
        prob_df = pd.DataFrame([[0 for _ in range(3)] for _ in range(4)],
                               index=['A', 'T', 'G', 'C'],
                               columns=[1, 2, 3])
        k, k_len = 0, 3
        for _ in range(self.tot_kmer()):
            _seq = self.seq[k:k + k_len]
            prob_df[1][_seq[0]] += 1
            prob_df[2][_seq[1]] += 1
            prob_df[3][_seq[2]] += 1
            k += k_len
        prob_df /= self.tot_kmer()
        # Gaussian correction
        return prob_df


if __name__ == '__main__':
    _seq = 'ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCAC'
    # _seq = 'ATGCGA'
    seq = Sequence(_seq, 11)
    print(seq.position_prob())
