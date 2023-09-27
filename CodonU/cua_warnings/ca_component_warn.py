from CodonU.cua_warnings.codon_usage_warns import CodonUsageWarning
import warnings


class CAComponentWarn(CodonUsageWarning):
    """
    Occurs when CA component is greater than data dimension
    """

    def __init__(self, comp, dim):
        self.msg = f'Analysis of {comp} components is not possible for your data with {dim} dimension. ' \
                   f'Analysing the data to maximum possible components.'

    def warn(self):
        warnings.warn(self.msg)
