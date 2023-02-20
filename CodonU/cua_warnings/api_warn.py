from .codon_usage_warns import CodonUsageWarning
import warnings


class ApiWarning(CodonUsageWarning):
    """
    Occurs when no API key is provided
    """

    def __init__(self):
        self.msg = 'No API key provided. Providing NCBI API key will increase the speed'
        super().__init__(self.msg)

    def warn(self):
        warnings.warn(self.msg)
