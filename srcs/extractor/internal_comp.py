from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature


def extract_cds_seq(seq: Seq, feature_location: FeatureLocation) -> Seq:
    """
    Extracts the CDS from a given sequence

    :param seq: Sequence
    :param feature_location: Feature location
    :return: The extracted feature
    """
    return feature_location.extract(seq)


def extract_prot_seq(feature: SeqFeature) -> Seq:
    """
    Returns the protein sequence reported in the report for the provided cds

    :param feature: The CDS
    :return: The protein sequence
    """
    return Seq(feature.qualifiers['translation'][0])
