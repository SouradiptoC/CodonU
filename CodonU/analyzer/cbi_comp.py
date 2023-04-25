from CodonU.analyzer.internal_comp import filter_reference, cbi
from warnings import filterwarnings
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.SeqIO import parse
import pandas as pd
from os.path import join, abspath
from CodonU.file_handler.internal_comp import is_file_writeable
from CodonU.file_handler import make_dir


def calculate_cbi(handle: str, genetic_code_num: int, min_len_threshold: int = 66, gene_analysis: bool = False,
                  save_file: bool = False, file_name: str = 'CBI_report', folder_path: str = 'Report') -> \
        dict[str, tuple[float, str] | dict[str, tuple[float, str]]]:
    """
    Calculates cbi values for each amino acid

    :param handle: Handle to the file, or the filename as a string
    :param genetic_code_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of nucleotide sequence to be considered as gene
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    :param save_file: Option for saving the values in xlsx format (Optional)
    :param file_name: Intended file name (Optional)
    :param folder_path: Folder path where image should be saved (optional)
    :return: The dictionary containing amino acid and cbi value, optimal codon pairs if gene_analysis is false,
    otherwise returns the dictionary containing gene name and dictionary containing amino acid and cbi value,
    optimal codon pairs
     """
    records = parse(handle, 'fasta')
    reference = filter_reference(records, min_len_threshold)
    filterwarnings('ignore')
    cbi_dict = dict()
    if gene_analysis:
        for i, seq in enumerate(reference):
            cbi_val_dict = dict()
            for aa in unambiguous_dna_by_id[genetic_code_num].protein_alphabet:
                cbi_val = cbi(aa, reference=[seq], genetic_code=genetic_code_num)
                cbi_val_dict.update({aa: cbi_val})
            cbi_dict.update({f'gene_{i + 1}': cbi_val_dict})
        if save_file:
            name = file_name + '.xlsx'
            make_dir(folder_path)
            file_path = join(folder_path, name)
            if is_file_writeable(file_path):
                df = pd.DataFrame.from_records(
                    [
                        (gene, aa, cbi_val[0], cbi_val[1])
                        for gene, cbi_vals in cbi_dict.items()
                        for aa, cbi_val in cbi_vals.items()
                    ],
                    columns=['Gene', 'AA', 'CBI_vals', 'Preferred_Codon']
                )
                df.to_excel(file_path, float_format='%.4f', columns=df.columns)
            print(f'The CBI score file can be found at: {abspath(file_path)}')
    else:
        for aa in unambiguous_dna_by_id[genetic_code_num].protein_alphabet:
            cbi_val = cbi(aa, reference, genetic_code_num)
            cbi_dict.update({aa: cbi_val})
        if save_file:
            name = file_name + '.xlsx'
            make_dir(folder_path)
            file_path = join(folder_path, name)
            if is_file_writeable(file_path):
                df = pd.DataFrame.from_records(
                    [
                        (aa, cbi_vals[0], cbi_vals[1])
                        for aa, cbi_vals in cbi_dict.items()
                    ],
                    columns=['AA', 'CBI_vals', 'Preferred_Codon']
                )
                df.to_excel(file_path, float_format='%.4f', columns=df.columns)
            print(f'The CBI score file can be found at: {abspath(file_path)}')
    return cbi_dict
