from .codon_usage_err import CodonUsageError


class NoProteinError(CodonUsageError):
    """
    Occurs when a complete category of amino acid based on sf values is not translated by the provided sequence
    """

    def __init__(self, seq):
        self.msg = 'A complete category of amino acid based on sf values is not translated by the provided sequence. ' \
                   'Try deleting the sequence or increase the threshold value for being considered as a gene\n. ' \
                   "For more details see 'The effective number of codons used in a gene' (1989). " \
                   f'The sequence is\n{seq}'
        super().__init__(self.msg)
