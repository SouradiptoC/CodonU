from typing import Tuple, Any

from Bio.SeqIO import parse
from Bio.Seq import Seq, MutableSeq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def extract_cds_lst(record: SeqRecord) -> tuple[Any, ...]:
    """
    Extracts the list of features if their type if CDS
    :param record: Original Sequence Record object from where the CDS is to be extracted
    :return: A tuple of FeatureLocation objects
    """
    cds_lst = [cds for cds in record.features if cds.type == 'CDS' and 'pseudo' not in cds.qualifiers.keys()]
    return tuple(cds_lst)


def extract_cds_seq(seq: Seq, feature_location: FeatureLocation) -> Seq:
    """
    Extracts the CDS from a given sequence
    :param seq: Sequence
    :param feature_location: Feature location
    :return: The extracted feature
    """
    return feature_location.extract(seq)


def extract_cds(record: SeqRecord, feature_location: FeatureLocation, cds_no: int = 0) -> SeqRecord:
    """
    Returns the CDS as a Sequence Record object
    :param record: Original Sequence Record object from where the CDS is to be extracted
    :param feature_location: The location of CDS
    :param cds_no: Number of CDS
    :return: The new Sequence Record object containing the CDS
    """
    cds = SeqRecord(
        seq=extract_cds_seq(record.seq, feature_location),
        id=f"{record.id} {record.annotations['organism']}",
        name=f"{record.annotations['organism']}",
        description=f"CDS_{cds_no}"
    )
    return cds


def extract_prot_seq(feature: SeqFeature) -> Seq:
    """
    Returns the protein sequence reported in the report for the provided cds
    :param feature: The CDS
    :return: The protein sequence
    """
    return Seq(feature.qualifiers['translation'][0])


def extract_prot(feature: SeqFeature, organism_name: str, cds_no: int = 0) -> SeqRecord:
    """
    Extracts protein sequences and return them for writing
    :param feature: The CDS
    :param organism_name: Name of the organism
    :param cds_no: Number of the CDS
    :return: The protein sequence suitable for being written is fasta format
    """
    if 'product' in feature.qualifiers.keys():
        description = f"{feature.qualifiers['product'][0]} CDS_{cds_no}"
    else:
        description = f"{feature.qualifiers['note'][0]} CDS_{cds_no}"
    prot = SeqRecord(
        seq=extract_prot_seq(feature),
        id=f"{feature.qualifiers['protein_id'][0]}",
        name=f"{organism_name}",
        description=description
    )
    return prot


def extract_exome(nuc_file_path: str, organism_name: str) -> SeqRecord:
    """
    Extracts the exome from given nucleotides
    :param nuc_file_path: The path to the nucleotide file
    :param organism_name: Name of the organism
    :return: The exome
    """
    records = parse(nuc_file_path, 'fasta')
    m_seq = MutableSeq('')
    for record in records:
        m_seq += record.seq[:-3]
    records = parse(nuc_file_path, 'fasta')
    *_, lst_record = records
    m_seq += lst_record.seq[-3:]
    exome = SeqRecord(
        seq=m_seq,
        id=lst_record.id,
        name=organism_name,
        description=f'whole exome of {organism_name}'
    )
    return exome
