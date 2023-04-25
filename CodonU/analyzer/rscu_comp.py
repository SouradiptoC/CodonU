from cai2 import RSCU
from CodonU.analyzer.internal_comp import filter_reference
from Bio.SeqIO import parse
from os.path import join, abspath
from CodonU.file_handler.internal_comp import is_file_writeable
from CodonU.file_handler import make_dir
import pandas as pd


def calculate_rscu(handle: str, genetic_code_num: int, min_len_threshold: int = 200, gene_analysis: bool = False,
                   save_file: bool = False, file_name: str = 'RSCU_report', folder_path: str = 'Report') -> \
        dict[str, float | dict[str, float]]:
    """
    Calculates rscu values for each codon

    :param handle: Handle to the file, or the filename as a string
    :param genetic_code_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    :param save_file: Option for saving the values in xlsx format (Optional)
    :param file_name: Intended file name (Optional)
    :param folder_path: Folder path where image should be saved (optional)
    :return: The dictionary containing codon and rscu value pairs if gene_analysis is false, otherwise the dictionary containing the gene name and the codon & rscu value pairs
    """
    records = parse(handle, 'fasta')
    references = filter_reference(records, min_len_threshold)
    if gene_analysis:
        rscu_dict = dict()
        for i, seq in enumerate(references):
            rscu_dict.update({f'gene_{i + 1}': RSCU([seq], genetic_code_num)})
        if save_file:
            name = file_name + '.xlsx'
            make_dir(folder_path)
            file_path = join(folder_path, name)
            if is_file_writeable(file_path):
                df = pd.DataFrame.from_records(
                    [
                        (gene, codon, rscu_val)
                        for gene, rscu_vals in rscu_dict.items()
                        for codon, rscu_val in rscu_vals.items()
                    ],
                    columns=['Gene', 'Codon', 'CAI_vals']
                )
                df.to_excel(file_path, float_format='%.4f', columns=df.columns)
            print(f'The RSCU score file can be found at: {abspath(file_path)}')
    else:
        reference = filter_reference(records, min_len_threshold)
        rscu_dict = RSCU(reference, genetic_code_num)
        if save_file:
            name = file_name + '.xlsx'
            make_dir(folder_path)
            file_path = join(folder_path, name)
            if is_file_writeable(file_path):
                df = pd.DataFrame.from_records(
                    [
                        (codon, rscu_val)
                        for codon, rscu_val in rscu_dict.items()
                    ],
                    columns=['Codon', 'CAI_vals']
                )
                df.to_excel(file_path, float_format='%.4f', columns=df.columns)
            print(f'The RSCU score file can be found at: {abspath(file_path)}')
    return rscu_dict
