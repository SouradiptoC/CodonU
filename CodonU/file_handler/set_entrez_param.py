from .internal_comp import set_entrez_email, set_entrez_api_key
from typing import Optional


def set_entrez_param(email: Optional[str] = None, api_key: Optional[str] = None) -> None:
    """
    Sets entrez parameters

    :param email: Email of the user (optional)
    :param api_key: API key of the user (optional)
    :raises EmailWarning: If no email is provided
    :raises ApiWarning: If no API key is provided
    """
    set_entrez_email(email)
    set_entrez_api_key(api_key)
