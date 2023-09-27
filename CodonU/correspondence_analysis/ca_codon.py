import pandas as pd
from CodonU.cua_warnings import CAComponentWarn
from CodonU.file_handler import make_dir
from CodonU.file_handler.internal_comp import is_file_writeable
from prince import CA
from typing import Optional
from os.path import join


def ca_codon(contingency_table: pd.DataFrame, n_components: int = 59,
             single_syn_codons: Optional[list[str]] = None, save_file: bool = False,
             file_name: str = 'CA_codon', folder_path: str = 'Report') -> CA:
    """
    Performs CA on codon frequency or RSCU contingency table\n
    **Note** Gene descriptions must be index headings and codons must column headings of the table

    :param contingency_table: The contingency table (pandas DataFrame object)
    :param n_components: Components for CA
    :param single_syn_codons: List of codons belonging to SF1
    :param save_file: Option for saving the values in xlsx format (Optional)
    :param file_name: Intended file name (Optional)
    :param folder_path: Folder path where image should be saved (optional)
    :return: The CA object
    """
    if single_syn_codons is None:
        # In standard codon table, ATG [-> Met] and TGG [-> Trp] belong to single synonymous codon family (SF1)
        single_syn_codons = ['ATG', 'TGG']

    # codons of SF1 needed to be dropped for CA
    _contingency_table = contingency_table.drop(single_syn_codons, axis=1)
    _contingency_table.replace(0, 0.000001, inplace=True)

    if n_components > len(_contingency_table.columns):
        warn = CAComponentWarn(n_components, len(_contingency_table.columns))
        warn.warn()
    ca = CA(random_state=42, n_components=n_components)
    ca.fit(_contingency_table)
    print(ca.eigenvalues_summary)

    if save_file:
        name = file_name + '.txt'
        make_dir(folder_path)
        file_path = join(folder_path, name)
        if is_file_writeable(file_path):
            with open(file_path, 'w') as out_file:
                print(ca.eigenvalues_summary, file=out_file)

    return ca
