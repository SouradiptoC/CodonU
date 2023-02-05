from .codon_usage_err import CodonUsageError


class FileNotEmptyError(CodonUsageError):
    """
    Occurs when a given file to write is not empty
    """

    def __init__(self, file_name=None):
        self.msg = f"{file_name} is not empty. Choose an empty file."
        super().__init__(self.msg)
