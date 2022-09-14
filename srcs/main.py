from Bio.Seq import Seq
import Bio
from Bio.Entrez import efetch
import Bio.SeqIO

Bio.Entrez.email = 'sourochaudhuri@gmail.com'

# strand = Seq('AATCG')
# print(strand.complement())
# i = len('MKQKEKEKIIKKTIVGGILGATMGYLSTSKSNPKSGTMTKDRLDHMKSLTSQLLNQQSQSKKSKSSVFSSFTKNKATKRKGSLWDKFKKDHKGDRKDKKSMKKTKINSRIFHNDDTSTKDSDNKGNVLNSLKKDKKEKTQSASKTPRTDAKKAKAKVKEEGSNHEKKNKGKTIFKKSEKSEKNGDQTKNEKARKKDAKKAEKNGKEAKITKEKEKPSLLDRFKRDDQSSKKTNKKKEKSNKDDSLKDSILNKFGKNGANKKKDLKKKKKKQKEKKGLLGKMRK')
# print(i)
# prot_1 = Seq('MKGKMRK')
# mRNA = prot_1.back_translate()
# len_mrna = len(mRNA)
# print(len_mrna)
k_len = 3
gc_str = 'gcGC'

codons_lst = list()
third_nucl = list()


def gc_count(seq: Bio.Seq.Seq, k_len: int = k_len):
    count = 0
    k = 0
    for _ in range(len(seq) // k_len):
        _seq = seq[k: k + k_len]
        # codons_lst.append(_seq)
        k += k_len
        if _seq[2] in gc_str:
            # third_nucl.append(_seq[2])
            count += 1
    return count


def gc_count_percentage(seq: Bio.Seq.Seq, k_len: int = k_len):
    count = gc_count(seq)
    tot_k_mer = len(seq) // k_len
    return (count / tot_k_mer) * 100


seq_1 = Bio.Seq.Seq('ATGAAGGGGAAGATGCGGAAG')
c = gc_count(seq_1)
# c_perc = gc_count_percentage(seq_1)
# print(c, c_perc)
# print(codons_lst)
# print(third_nucl)
#
# codons_lst.clear()
# third_nucl.clear()

handle = efetch(db='nucleotide', id='AF3335109.1', rettype='fasta')
record = Bio.SeqIO.read(handle, 'fasta')
print(len(record))
c_2 = gc_count(record.seq)
c_2_perc = gc_count_percentage(record.seq)
# print(codons_lst)
# print(third_nucl)
print(c_2, c_2_perc)
