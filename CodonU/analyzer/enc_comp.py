from warnings import filterwarnings
from CodonU.analyzer.internal_comp import filter_reference, enc
from Bio.SeqIO import parse
import pandas as pd
from os.path import join, abspath
from CodonU.file_handler.internal_comp import is_file_writeable
from CodonU.file_handler import make_dir


def calculate_enc(handle: str, genetic_code_num: int, min_len_threshold=200, gene_analysis: bool = False,
                  save_file: bool = False, file_name: str = 'ENc_report', folder_path: str = 'Report') -> \
        float or dict[str, float]:
    """
    Calculates ENc value for a given sequences

    :param handle: Handle to the file, or the filename as a string
    :param genetic_code_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    :param save_file: Option for saving the values in xlsx format (Optional)
    :param file_name: Intended file name (Optional)
    :param folder_path: Folder path where image should be saved (optional)
    :return: The ENc value if gene_analysis is false, else a dictionary containing gene number and corresponding ENc value
    """
    records = parse(handle, 'fasta')
    filtered_records = filter_reference(records, min_len_threshold, 'nuc')
    filterwarnings('ignore')
    if gene_analysis:
        enc_dict = dict()
        for record in filtered_records:
            enc_dict.update({record.description: enc([record.seq], genetic_code_num)})
        if save_file:
            name = file_name + '.xlsx'
            make_dir(folder_path)
            file_path = join(folder_path, name)
            if is_file_writeable(file_path):
                df = pd.DataFrame(enc_dict.items(), columns=['Gene', 'ENc_val'])
                df.to_excel(file_path, float_format='%.4f', columns=df.columns)
            print(f'The ENc score file can be found at: {abspath(file_path)}')
        return enc_dict
    else:
        sequences = [record.seq for record in filtered_records]
        if save_file:
            name = file_name + '.xlsx'
            make_dir(folder_path)
            file_path = join(folder_path, name)
            if is_file_writeable(file_path):
                df = pd.DataFrame({'Whole Genome': enc(sequences, genetic_code_num)}.items(),
                                  columns=['Genome', 'ENc_vals'])
                df.to_excel(file_path, float_format='%.4f', columns=df.columns)
            print(f'The ENc score file can be found at: {abspath(file_path)}')
        return enc(sequences, genetic_code_num)
