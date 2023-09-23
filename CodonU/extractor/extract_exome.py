from Bio.Seq import MutableSeq
from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord


def extract_exome(nuc_file_path: str, organism_name: str, exclude_stops: bool = True) -> SeqRecord:
    """
    Extracts the exome from given nucleotides

    :param nuc_file_path: The path to the nucleotide file
    :param organism_name: Name of the organism
    :param exclude_stops: If true, intermediate stops codons are excluded from exome
    :return: The exome
    """
    records = parse(nuc_file_path, 'fasta')
    m_seq = MutableSeq('')
    if exclude_stops:
        for record in records:
            m_seq += record.seq[:-3]
        records = parse(nuc_file_path, 'fasta')
        *_, lst_record = records
        m_seq += lst_record.seq[-3:]
    else:
        for record in records:
            m_seq += record.seq
        records = parse(nuc_file_path, 'fasta')
        *_, lst_record = records
    exome = SeqRecord(
        seq=m_seq,
        id=lst_record.id,
        name=organism_name,
        description=f'whole exome of {organism_name}'
    )
    return exome
