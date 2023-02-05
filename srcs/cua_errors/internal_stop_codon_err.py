from .codon_usage_err import CodonUsageError


class InternalStopCodonError(CodonUsageError):
    """
    Occurs when internal stop codon is present in the genome
    """

    def __init__(self, position):
        """
        :param position: Position of detected nucleotide
        """
        self.msg = f'Internal stop codon detected at position {position + 1}. Enter valid genome.'
        super().__init__(self.msg)
