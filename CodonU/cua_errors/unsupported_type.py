from .codon_usage_err import CodonUsageError


class UnsupportedType(CodonUsageError):
    def __init__(self):
        self.msg = 'The type of sequence is not supported. ' \
                   "Supported types are 'nucleotide' and 'aa'.\n" \
                   'For raising the issue about the compatibility of other types of sequence analysis, ' \
                   'please visit https://github.com/SouradiptoC/CodonU/issues'
        super().__init__(self.msg)
