from .codon_usage_warns import CodonUsageWarning
import warnings


class MissingCodonWarning(CodonUsageWarning):
    """
    Occurs when no codon in the given reference sequence list translates to a certain amino acid
    """

    def __init__(self, aa: str):
        self.msg = f'No codon in the given reference sequence list translates to {aa}'
        super().__init__(self.msg)

    def warn(self):
        warnings.warn(self.msg)
