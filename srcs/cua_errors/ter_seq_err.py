from .codon_usage_err import CodonUsageError


class TerSeqError(CodonUsageError):
    """
    Occurs when terminal codons are not valid
    """

    def __init__(self, code):
        """
        :param code: 1 for starting, -1 for ending
        """
        msg_dict = {
            1: 'The first codon is not a valid starting codon.',
            -1: 'The last codon is not a valid ending codon.'
        }
        self.msg = msg_dict[code]
        super().__init__(self.msg)
