from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from CodonU.extractor.internal_comp import extract_prot_seq


def extract_prot(feature: SeqFeature, organism_name: str) -> SeqRecord:
    """
    Extracts protein sequences and return them for writing

    :param feature: The CDS
    :param organism_name: Name of the organism
    :param cds_no: Number of the CDS
    :return: The protein sequence suitable for being written is fasta format
    """
    if 'product' in feature.qualifiers.keys():
        description = f"{feature.qualifiers['product'][0]}"
    else:
        description = f"{feature.qualifiers['note'][0]}"
    prot = SeqRecord(
        seq=extract_prot_seq(feature),
        id=f"{feature.qualifiers['protein_id'][0]}",
        name=f"{organism_name}",
        description=description
    )
    return prot
