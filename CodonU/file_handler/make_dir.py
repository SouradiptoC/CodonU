import os
from CodonU.cua_logger import *


def make_dir(path: str) -> None:
    """
    Makes a directory if not present already

    :param path: Path of the directory
    """
    if not os.path.isdir(path):
        os.mkdir(path)
        name = path.split('/')
        console_log.info(f"{name[-1]} created successfully")
        file_log.info(f"{name[-1]} created successfully")
