class CodonUsageWarning(Warning):
    """
    Accounts for all warnings while executing the CodonU
    """

    def __init__(self, msg=None):
        self.msg = msg
