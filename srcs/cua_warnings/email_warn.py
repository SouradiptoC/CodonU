from .codon_usage_warns import CodonUsageWarning
import warnings


class EmailWarning(CodonUsageWarning):
    """
    Occurs when no email id is provided
    """

    def __init__(self):
        self.msg = 'No Email provided. Providing email will increase the speed'
        super().__init__(self.msg)

    def warn(self):
        warnings.warn(self.msg)
