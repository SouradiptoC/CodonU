from Bio import Entrez
from Bio.SeqIO import read
from Bio.SeqRecord import SeqRecord


def get_gb(accession_id: str) -> SeqRecord:
    """
    Gets the Sequence Record object from a given accession number

    :param accession_id: Provided accession number
    :return: The Sequence Record object
    """
    print(f"Retrieval started")
    handle = Entrez.efetch(db='nucleotide', id=accession_id, rettype='gb', retmode='text')
    record = read(handle, 'gb')
    print(f"Genbank file of {record.annotations['source']} retrieved successfully")
    return record
