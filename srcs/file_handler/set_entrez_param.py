from .internal_comp import set_entrez_email, set_entrez_api_key


def set_entrez_param(email: str | None = None, api_key: str | None = None) -> None:
    """
    Sets entrez parameters

    :param email: Email of the user (optional)
    :param api_key: API key of the user (optional)
    :raises EmailWarning: If no email is provided
    :raises ApiWarning: If no API key is provided
    """
    set_entrez_email(email)
    set_entrez_api_key(api_key)
