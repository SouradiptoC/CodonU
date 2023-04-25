from cai2 import CAI
from warnings import filterwarnings
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.SeqIO import parse
from CodonU.analyzer.internal_comp import filter_reference
import pandas as pd
from os.path import join, abspath
from CodonU.file_handler.internal_comp import is_file_writeable
from CodonU.file_handler import make_dir


def calculate_cai(handle: str, genetic_code_num: int, min_len_threshold: int = 200, gene_analysis: bool = False,
                  save_file: bool = False, file_name: str = 'CAI_report', folder_path: str = 'Report') -> \
        dict[str, float | dict[str, float]]:
    """
    Calculates cai values for each codon

    :param handle: Handle to the file, or the filename as a string
    :param genetic_code_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    :param save_file: Option for saving the values in xlsx format (Optional)
    :param file_name: Intended file name (Optional)
    :param folder_path: Folder path where image should be saved (optional)
    :return: The dictionary containing codon and cai value pairs if gene_analysis is False, otherwise returns the
    dictionary containing gene name and corresponding codon and cai value pairs
    """
    filterwarnings('ignore')
    cai_dict = dict()
    records = parse(handle, 'fasta')
    reference = filter_reference(records, min_len_threshold)
    if gene_analysis:
        for i, seq in enumerate(reference):
            cai_val_dict = dict()
            for codon in unambiguous_dna_by_id[genetic_code_num].forward_table:
                cai_val = CAI(codon, reference=[seq], genetic_code=genetic_code_num)
                cai_val_dict.update({codon: cai_val})
            cai_dict.update({f'gene_{i + 1}': cai_val_dict})
        if save_file:
            name = file_name + '.xlsx'
            make_dir(folder_path)
            file_path = join(folder_path, name)
            if is_file_writeable(file_path):
                df = pd.DataFrame.from_records(
                    [
                        (gene, codon, cai_val)
                        for gene, cai_vals in cai_dict.items()
                        for codon, cai_val in cai_vals.items()
                    ],
                    columns=['Gene', 'Codon', 'CAI_vals']
                )
                df.to_excel(file_path, float_format='%.4f', columns=df.columns)
            print(f'The CAI score file can be found at: {abspath(file_path)}')
    else:
        for codon in unambiguous_dna_by_id[genetic_code_num].forward_table:
            cai_val = CAI(codon, reference=reference, genetic_code=genetic_code_num)
            cai_dict.update({codon: cai_val})
        if save_file:
            name = file_name + '.xlsx'
            make_dir(folder_path)
            file_path = join(folder_path, name)
            if is_file_writeable(file_path):
                df = pd.DataFrame(cai_dict.items(), columns=['Codon', 'CAI_vals'])
                df.to_excel(file_path, float_format='%.4f', columns=df.columns)
            print(f'The CAI score file can be found at: {abspath(file_path)}')
    return cai_dict
