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
