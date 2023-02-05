from typing import Any

from Bio.SeqRecord import SeqRecord


def extract_cds_lst(record: SeqRecord) -> tuple[Any, ...]:
    """
    Extracts the list of features if their type is CDS

    :param record: Original Sequence Record object from where the CDS is to be extracted
    :return: A tuple of FeatureLocation objects
    """
    cds_lst = [cds for cds in record.features if cds.type == 'CDS' and 'pseudo' not in cds.qualifiers.keys()]
    return tuple(cds_lst)
