from Bio.SeqIO import parse

from math import isnan
from Bio.Seq import Seq
from Bio.Data.IUPACData import unambiguous_dna_letters

from Bio.SeqUtils import GC, GC123
from CAI import CAI

from warnings import filterwarnings


def gc_123(seq: Seq | str) -> tuple:
    return GC123(seq)


def cai(reference: list, genetic_code: int):
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
    # print(gc_123('GGG'))
    handle = '/home/souro/Projects/final_yr/Results/Nucleotide/Staphylococcus_agnetis_nucleotide.fasta'
    lst = []
    records = parse(handle, 'fasta')
    for record in records:
        lst.append(record.seq)
    cai_dct = cai(lst, 11)
    print(cai_dct.keys())
