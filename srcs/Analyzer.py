from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def extract_cds_lst(record: SeqRecord) -> tuple:
    """
    Extracts the list of features if their type if CDS
    :param record: Original Sequence Record object from where the CDS is to be extracted
    :return: A tuple of FeatureLocation objects
    """
    # TODO implement len of gene >= 300 or find 25% of median
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
    return Seq(feature.qualifiers['translation'][0])


def extract_prot(feature: SeqFeature, organism_name: str, cds_no: int = 0) -> SeqRecord:
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


if __name__ == '__main__':
    pass
