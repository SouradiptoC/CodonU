from .codon_usage_warns import CodonUsageWarning
import warnings


class NoCodonWarning(CodonUsageWarning):
    """
    Occurs when a codon is not present in a given sequence
    """

    def __init__(self, codon: str):
        """
        :param codon: The codon which is absent
        """
        self.msg = f'{codon} is not present in given sequence'
        super().__init__(self.msg)

    def warn(self):
        warnings.warn(self.msg)
