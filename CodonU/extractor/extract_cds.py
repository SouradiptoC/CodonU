from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord

from .internal_comp import extract_cds_seq


def extract_cds(record: SeqRecord, feature_location: FeatureLocation) -> SeqRecord:
    """
    Returns the CDS as a Sequence Record object

    :param record: Original Sequence Record object from where the CDS is to be extracted
    :param feature_location: The location of CDS
    :return: The new Sequence Record object containing the CDS
    """
    cds = SeqRecord(
        seq=extract_cds_seq(record.seq, feature_location),
        id=f"{record.id} {record.annotations['organism']}",
        name=f"{record.annotations['organism']}",
        description=f"{record.description}"
    )
    return cds
