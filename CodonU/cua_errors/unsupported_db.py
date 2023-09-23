from .codon_usage_err import CodonUsageError


class UnsupportedDatabase(CodonUsageError):
    def __init__(self):
        self.msg = 'The Database is not supported. ' \
                   "Supported databases are 'tRNADB_CE' and 'GtRNAdb'.\n" \
                   'For raising the issue about the compatibility of other databases, ' \
                   'please visit https://github.com/SouradiptoC/CodonU/issues'
        super().__init__(self.msg)
