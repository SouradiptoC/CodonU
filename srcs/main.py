import Bio
import Bio.SeqIO
from Bio.Seq import Seq
from Bio.Entrez import efetch

Bio.Entrez.email = 'sourochaudhuri@gmail.com'


def tot_k_mers(seq: Bio.Seq.Seq, k_len: int = 3):
    return len(seq) // k_len


def gc_count(seq: Bio.Seq.Seq, k_len: int = 3):
    gc_str = 'GCgc'
    count = 0
    k = 0
    for _ in range(tot_k_mers(seq)):
        _seq = seq[k: k + k_len]
        k += k_len
        if _seq[2] in gc_str:
            count += 1
    return count


def gc_count_percentage(seq: Bio.Seq.Seq):
    count = gc_count(seq)
    tot_k_mer = tot_k_mers(seq)
    return (count / tot_k_mer) * 100


if __name__ == '__main__':
    pass
