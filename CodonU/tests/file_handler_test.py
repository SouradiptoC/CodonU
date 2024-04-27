from CodonU.file_handler import *
from Bio import Entrez
from Bio import SeqIO


def test_get_gb():
    res = get_gb('AP009351')
    handle = Entrez.efetch(db='nucleotide', id='AP009351', rettype='gb', retmode='text')
    record = SeqIO.read(handle, 'gb')
    assert res.id == record.id, "get_gb is not working"
