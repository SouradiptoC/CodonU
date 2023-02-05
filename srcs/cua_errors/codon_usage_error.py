class CodonUsageError(Exception):
    """
    Accounts for all error that can happen during executing the files
    """

    def __init__(self, msg=None):
        self.msg = msg

    def __str__(self):
        return self.msg
