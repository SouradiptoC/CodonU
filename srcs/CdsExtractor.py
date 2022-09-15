from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord


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


if __name__ == '__main__':
    pass
