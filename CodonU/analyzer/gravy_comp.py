from warnings import filterwarnings
from CodonU.analyzer.internal_comp import filter_reference, gravy
from Bio.SeqIO import parse
from CodonU.file_handler.internal_comp import is_file_writeable
from CodonU.file_handler import make_dir
import pandas as pd
from os.path import join, abspath


def calculate_gravy(handle: str, min_len_threshold: int = 66, gene_analysis: bool = False, save_file: bool = False,
                    file_name: str = 'GRAVY_report', folder_path: str = 'Report') -> dict[str, float] | float:
    """
    Calculates the gravy score for a given protein sequence

    :param handle: Handle to the file, or the filename as a string
    :param min_len_threshold: Minimum length of protein sequence to be considered as gene
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    :param save_file: Option for saving the values in xlsx format (Optional)
    :param file_name: Intended file name (Optional)
    :param folder_path: Folder path where image should be saved (optional)
    :return: The GRAVY score of given sequence if gene_analysis is false, else the dictionary containing gene number and
    corresponding GRAVY score
    """
    records = parse(handle, 'fasta')
    references = filter_reference(records, min_len_threshold)
    filterwarnings('ignore')
    if gene_analysis:
        gravy_dict = dict()
        for i, seq in enumerate(references):
            gravy_dict.update({f'prot_seq{i + 1}': gravy(seq)})
        if save_file:
            name = file_name + '.xlsx'
            make_dir(folder_path)
            file_path = join(folder_path, name)
            if is_file_writeable(file_path):
                df = pd.DataFrame(gravy_dict.items(), columns=['Protein_name', 'Gravy_score'])
                df.to_excel(file_path, float_format='%.4f', columns=df.columns)
            print(f'The GRAVY score file can be found at: {abspath(file_path)}')
        return gravy_dict
    else:
        seq = ''.join([str(_seq) for _seq in references])
        if save_file:
            name = file_name + '.xlsx'
            make_dir(folder_path)
            file_path = join(folder_path, name)
            if is_file_writeable(file_path):
                df = pd.DataFrame({'Prot_seq': gravy(seq)}.items(), columns=['Protein_name', 'Gravy_score'])
                df.to_excel(file_path, float_format='%.4f', columns=df.columns)
            print(f'The GRAVY score file can be found at: {abspath(file_path)}')
        return gravy(seq)
