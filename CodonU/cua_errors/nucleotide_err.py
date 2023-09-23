from .codon_usage_err import CodonUsageError


class NucleotideError(CodonUsageError):
    """
    Occurs when an ambiguous or invalid nucleotide is present in genome
    """

    def __init__(self, code):
        """
        :param code: 1 for ambiguous, 2 for invalid
        """
        if code == 1:
            self.msg = f"Ambiguous nucleotide detected. Provide an unambiguous genome."
        elif code == 2:
            self.msg = f"Invalid nucleotide detected. Provide a valid genome."
