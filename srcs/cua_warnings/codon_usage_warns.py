class CodonUsageWarning(Warning):
    """
    Accounts for all warnings while executing the srcs
    """

    def __init__(self, msg=None):
        self.msg = msg
