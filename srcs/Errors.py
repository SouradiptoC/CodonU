import warnings


class CodonUsageError(Exception):
    """
    Accounts for all error that can happen during executing the files
    """

    def __init__(self, msg=None):
        self.msg = msg

    def __str__(self):
        return self.msg


class FileNotEmptyError(CodonUsageError):
    """
    Occurs when a given file to write is not empty
    """

    def __init__(self, file_name=None):
        self.msg = f"{file_name} is not empty. Choose an empty file."
        super().__init__(self.msg)


class InternalStopCodonError(CodonUsageError):
    """
    Occurs when internal stop codon is present in the genome
    """

    def __init__(self, position):
        """
        :param position: Position of detected nucleotide
        """
        self.msg = f'Internal stop codon detected at position {position + 1}. Enter valid genome.'
        super().__init__(self.msg)


class NucleotideError(CodonUsageError):
    """
    Occurs when an ambiguous or invalid nucleotide is present in genome
    """

    def __init__(self, code, position):
        """
        :param code: 1 for ambiguous, 2 for invalid
        :param position: position of detected nucleotide
        """
        if code == 1:
            self.msg = f"Ambiguous nucleotide detected at position {position + 1}. Provide an unambiguous genome."
        elif code == 2:
            self.msg = f"Invalid nucleotide detected at position {position + 1}. Provide a valid genome."


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


class ThresholdError(CodonUsageError):
    """
    Occurs when the provided threshold is not in limit
    """

    def __init__(self, seq):
        self.msg = 'A complete category of amino acid based on sf values is not translated by the provided sequence. ' \
                   'Try deleting the sequence or increase the threshold value for being considered as a gene\n. ' \
                   "For more details see 'The effective number of codons used in a gene' (1989). " \
                   f'The sequence is\n{seq}'
        super().__init__(self.msg)


##############################################################################################
class CodonUsageWarning(Exception):
    """
    Accounts for all warnings while executing the srcs
    """

    def __init__(self, msg=None):
        self.msg = msg


class EmailWarning(CodonUsageWarning):
    """
    Occurs when no email id is provided
    """

    def __init__(self):
        self.msg = 'No Email provided. Providing email will increase the speed'
        super().__init__(self.msg)

    def warn(self):
        warnings.warn(self.msg)


class ApiWarning(CodonUsageWarning):
    """
    Occurs when no API key is provided
    """

    def __init__(self):
        self.msg = 'No API key provided. Providing NCBI API key will increase the speed'
        super().__init__(self.msg)

    def warn(self):
        warnings.warn(self.msg)


class MissingCodonWarning(CodonUsageWarning):
    """
    Occurs when no codon in the given reference sequence list translates to a certain amino acid
    """

    def __init__(self, aa: str):
        self.msg = f'No codon in the given reference sequence list translates to {aa}'
        super().__init__(self.msg)

    def warn(self):
        warnings.warn(self.msg)


class NoSynonymousCodonWarning(CodonUsageWarning):
    """
    Occurs when only one codon in the given reference sequence list translates to a certain amino acid
    """

    def __init__(self, aa):
        self.msg = f'There is only one codon for {aa}'

    def warn(self):
        warnings.warn(self.msg)
