from .codon_usage_err import CodonUsageError


class BadSequenceError(CodonUsageError):
    """
    Occurs when the sequence is bad i.e. length of the sequence is not divisible by 3
    """

    def __init__(self, seq, code):
        """
        :param seq: The nucleotide sequence which is bad
        """
        if code == 1:
            self.msg = f"Length of the sequence is not divisible by 3. Check the sequence below.\n" \
                       f"{seq}"
        if code == 2:
            self.msg = f"Length of the sequence is not divisible by 3. Check the sequence below.\n" \
                       f"{seq.description}: {seq.seq}"
