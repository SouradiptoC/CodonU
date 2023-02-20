from .codon_usage_warns import CodonUsageWarning
import warnings


class NoSynonymousCodonWarning(CodonUsageWarning):
    """
    Occurs when only one codon in the given reference sequence list translates to a certain amino acid
    """

    def __init__(self, aa):
        self.msg = f'There is only one codon for {aa}'

    def warn(self):
        warnings.warn(self.msg)
