from .codon_usage_err import CodonUsageError


class NucleotideError(CodonUsageError):
    """
    Occurs when an ambiguous or invalid nucleotide is present in genome
    """

    def __init__(self, code, position):
        """
        :param code: 1 for ambiguous, 2 for invalid
        :param position: position of detected nucleotide
        """
        if code == 1:
            self.msg = f"Ambiguous nucleotide detected at position {position + 1}. Provide an unambiguous genome."
        elif code == 2:
            self.msg = f"Invalid nucleotide detected at position {position + 1}. Provide a valid genome."
