import os
from Bio import Entrez
from CodonU.cua_errors import FileNotEmptyError
from CodonU.cua_warnings import EmailWarning, ApiWarning


def set_entrez_email(email: str | None) -> None:
    """
    Sets Bio.Entrez.email parameter to given email

    :param email: Email of user
    :raises EmailWarning: If no email is provided
    """
    if email:
        print('Setting provided email to entrez.email')
        Entrez.email = email
    else:
        warning = EmailWarning()
        warning.warn()


def set_entrez_api_key(api_key: str | None) -> None:
    """
    Sets Bio.Entrez.api_key parameter to given api_key

    :param api_key: API key of the user
    :raises ApiWarning: If no API key is provided
    """
    if api_key:
        print('Setting provided API key to entrez.api_key')
        Entrez.api_key = api_key
    else:
        warning = ApiWarning()
        warning.warn()


def is_file_empty(path: str) -> bool:
    """
    Checks if an existing file is empty

    :param path: Path to the file
    :return: True if empty else false
    :raises FileNotEmptyError: If the given file to write is not empty
    """
    if os.stat(path).st_size == 0:
        return True
    else:
        name = path.split('/')
        dec = input(f"{name[-1]} already exists. Want to re-write (y/n): ")
        if dec == 'y':
            return True
        raise FileNotEmptyError(path)


def is_file(path: str) -> bool:
    """
    Checks if file exists or not

    :param path: Path to the file
    :return: True if exists else False
    """
    return os.path.isfile(path)


def is_file_writeable(path: str):
    if is_file(path) and os.stat(path).st_size != 0:
        flg = input(
            'Provided file not empty! Your action will result into completely changing the content of the file. Proceed [y/n]?: ')
        if flg in ['y', 'Y']:
            return True
        raise FileNotEmptyError(path)
    else:
        return True
