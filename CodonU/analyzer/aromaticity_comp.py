from CodonU.analyzer.internal_comp import filter_reference, aromaticity
from CodonU.file_handler.internal_comp import is_file_writeable
from CodonU.file_handler import make_dir
from Bio.SeqIO import parse
import pandas as pd
from os.path import join, abspath


def calculate_aromaticity(handle: str, min_len_threshold: int = 66, gene_analysis: bool = False,
                          save_file: bool = False, file_name: str = 'Aroma_report', folder_path: str = 'Report') -> \
        dict[str, float] | float:
    """
    Calculates the aromaticity score for a given protein sequence

    :param handle: Handle to the file, or the filename as a string
    :param min_len_threshold: Minimum length of protein sequence to be considered as gene
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    :param save_file: Option for saving the values in xlsx format (Optional)
    :param file_name: Intended file name (Optional)
    :param folder_path: Folder path where image should be saved (optional)
    :return: The aromaticity score of given sequence if gene_analysis is false, else the dictionary containing
    gene number and corresponding GRAVY score
    """
    records = parse(handle, 'fasta')
    references = filter_reference(records, min_len_threshold)
    if gene_analysis:
        aroma_dict = dict()
        for i, seq in enumerate(references):
            aroma_dict.update({f'prot_seq{i + 1}': aromaticity(seq)})
        if save_file:
            name = file_name + '.xlsx'
            make_dir(folder_path)
            file_path = join(folder_path, name)
            if is_file_writeable(file_path):
                df = pd.DataFrame(aroma_dict.items(), columns=['Protein_name', 'Aroma_score'])
                df.to_excel(file_path, float_format='%.4f', columns=df.columns)
            print(f'The Aromaticity score file can be found at: {abspath(file_path)}')
        return aroma_dict
    else:
        seq = ''.join([str(_seq) for _seq in references])
        if save_file:
            name = file_name + '.xlsx'
            make_dir(folder_path)
            file_path = join(folder_path, name)
            if is_file_writeable(file_path):
                df = pd.DataFrame({'Prot_seq': aromaticity(seq)}.items(), columns=['Protein_name', 'Aroma_score'])
                df.to_excel(file_path, float_format='%.4f', columns=df.columns)
            print(f'The Aromaticity score file can be found at: {abspath(file_path)}')
        return aromaticity(seq)
