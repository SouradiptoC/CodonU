from .codon_usage_err import CodonUsageError


class CodonTableExistsError(CodonUsageError):
    """
    Occurs when id, name or alt_name of a new table is same with existing tables
    """

    def __init__(self, code, val):
        """
        :param code: 1 for id, 2 for name, 3 for alt_name
        :param val: The id or name or alt_name provided by user
        """
        code_msg = {
            1: "id", 2: "name", 3: "alternative name"
        }
        self.msg = f"The {code_msg[code]} [{val}] already is in use for standard tables by NCBI. " \
                   f"Use the available table or give other values for the {code_msg[code]}.\n" \
                   f"Check standard tables at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"
