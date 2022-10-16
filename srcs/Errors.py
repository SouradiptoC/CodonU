class CodonUsageError(Exception):
    """
    Accounts for all error that can happen during executing the files
    """

    def __init__(self, msg=None):
        self.msg = msg

    def __str__(self):
        return self.msg


class FileNotEmpty(CodonUsageError):
    """
    Occurs when a given file to write is not empty
    """

    def __init__(self, file_name=None):
        self.file_name = file_name
        self.msg = f"{self.file_name} is not empty. Choose an empty file."
        super().__init__(self.msg)


class NoEmailError(CodonUsageError):
    """
    Occurs when no email is provided for Bio.Entrez.email parameter
    """

    def __init__(self):
        self.msg = f"No email provided for accessing NCBI utilities. Provide a valid email to access."
        super().__init__(self.msg)


class TerSeqError(CodonUsageError):
    """
    Occurs when terminal codons are not valid
    """

    def __init__(self, code):
        msg_dict = {
            1: 'The first codon is not a valid starting codon',
            -1: 'The last codon is not a valid ending codon'
        }
        self.msg = msg_dict[code]
        super().__init__(self.msg)
