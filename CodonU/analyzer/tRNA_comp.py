import pandas as pd
import warnings
import os
from typing import Optional
from Bio.SeqIO import parse
from gtAI.bygaft import gene_algo_corr
from gtAI.gtAI import tRNADB_CE, GtRNAdb, dict_codon_anticodon, dict_codon_anticodon_count, abs_Wi, rel_Wi, calc_Tai
from gtAI.new_flow import ENc_calc, ENc_calc_ref, RSCU_calc
from CodonU.analyzer.internal_comp import not_contains_amb_letter
from CodonU.file_handler import make_dir
from CodonU.file_handler.internal_comp import is_file, is_file_writeable
from CodonU.cua_errors import UnsupportedDatabase

# for ignoring future warnings
warnings.filterwarnings('ignore')


def get_anticodon_count_dict(url: str, database: str) -> dict[str, int]:
    """
    Retrieves the anticodon table from given link\n
    **NOTE:** The *database* can have only two values, i.e. "*tRNADB_CE*" and "*GtRNAdb*"
        * For using *tRNADB_CE*, please visit http://trna.ie.niigata-u.ac.jp/cgi-bin/trnadb/index.cgi
        * For using *GtRNAdb*, please visit http://gtrnadb.ucsc.edu/

    :param url: URL to anticodon table
    :param database: Type of database from the above options
    :return: The dictionary containing anticodon as key and count as val
    :raise UnsupportedDatabase: If *database* has other values than mentioned
    """
    if database == "tRNADB_CE":
        return tRNADB_CE(url)
    elif database == "GtRNAdb":
        return GtRNAdb(url)
    else:
        raise UnsupportedDatabase()


def calculate_gtai(handle: str, anticodon_dict: dict, genetic_code_num: int, reference: Optional[str] = None,
                   size_pop: int = 60, generation_num: int = 100,
                   save_file: bool = False, file_name: str = 'tAI_report', folder_path: str = 'Report') -> \
        tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Calculates the gtAI value for each gene according to Anwar et al., 2023\n
    The function returns following dataframes:
        * tai_df: The dataframe contains gene description and tAI values
        * abs_wi_df: The dataframe contains each anticodon and absolute weights according to the paper
        * rel_wi_df: The dataframe contains each anticodon and relative weights according to the paper
    \n**Note:** The function will generate a file named 'best_fit.py'

    :param handle: Path to the fasta file as a string
    :param anticodon_dict: The dictionary containing anticodon as key and count as value
    :param genetic_code_num: Genetic table number for codon table
    :param reference: Path to the reference fasta file as a string (Optional)
    :param size_pop: A parameter for the genetic algorithm to identify the population size (Optional)
    :param generation_num: A parameter for the genetic algorithm to identify the generation number (Optional)
    :param save_file: Option for saving the values in xlsx format (Optional)
    :param file_name: Intended file name (Optional)
    :param folder_path: Folder path where image should be saved (optional)
    :return: A tuple of 3 dataframes, as discussed earlier
    :raises FileExistsError: If re-write permission is not given for the file `best_fit.py`
    :raises ImportError: If `best_fit.py` is not created or deleted after creation
    """
    if reference is None:
        enc_low_values = ENc_calc(handle)
        rscu_values = RSCU_calc(handle, enc_low_values, genetic_code_num)
    else:
        enc_low_values = ENc_calc_ref(reference)
        rscu_values = RSCU_calc(reference, enc_low_values, genetic_code_num)

    bac_flg = True if genetic_code_num == 11 else False  # as genetic codon table 11 corresponds to bacteria table

    file_path = os.path.join(os.getcwd(), 'best_fit.py')
    if is_file(file_path):
        ans = input("'best_fit.py' already exists! Want to re-write(y/n): ")
        if ans in 'Yy':
            os.remove(file_path)
        else:
            raise FileExistsError("'best_fit.py' exists!")
    gene_algo_corr(anticodon_dict, genetic_code_num, rscu_values, bac_flg, size_pop, generation_num)
    # this command will create a file named best_fit.py
    # if best_fit.py exists, then it must be deleted and is done by the block of codes earlier

    while True:
        try:
            import best_fit
            best_fit_ = best_fit.best_fit
            _, variants, _ = map(list, zip(*best_fit_))
            Sij = variants[-1]
            Sug, Sci, Sai, Sgu, Sal = Sij[:5]
            break
        except ImportError as ie:
            raise ie
        except Exception as e:
            raise e

    codon_anticodon_dict = dict_codon_anticodon(anticodon_dict)
    codon_anticodon_count_dict = dict_codon_anticodon_count(codon_anticodon_dict, anticodon_dict, bacteria=bac_flg)
    abs_wi_dict = abs_Wi(codon_anticodon_count_dict, Sug=Sug, Sci=Sci, Sai=Sai, Sgu=Sgu, Sal=Sal, bacteria=bac_flg)
    rel_wi_dict = rel_Wi(abs_wi_dict, genetic_code_number=genetic_code_num)

    tai_dict = dict()
    records = parse(handle, 'fasta')
    for record in records:
        seq = str(record.seq).upper()
        if not_contains_amb_letter(seq):
            tai = calc_Tai(seq, rel_dict_wi=rel_wi_dict, genetic_code_number=genetic_code_num)
            tai_dict.update({record.description: tai})
        else:
            pass

    tai_df = pd.DataFrame(list(tai_dict.items()))
    tai_df.columns = ["gene_description", "tAI_values"]
    abs_wi_df = pd.DataFrame(list(abs_wi_dict.items()))
    abs_wi_df.columns = ["anti_codon", "abs_wi_value"]
    rel_wi_df = pd.DataFrame(list(rel_wi_dict.items()))
    rel_wi_df.columns = ["anti_codon", "rel_wi_value"]

    if save_file:
        names = ["", "_abs_wi", "_rel_wi"]
        dfs = [tai_df, abs_wi_df, rel_wi_df]
        for i in range(len(names)):
            name = file_name + names[i] + '.xlsx'
            make_dir(folder_path)
            file_path = os.path.join(folder_path, name)
            if is_file_writeable(file_path):
                dfs[i].to_excel(file_path, float_format='%.4f', columns=dfs[i].columns)
            print(f'The {name} file can be found at: {os.path.abspath(file_path)}')

    return tai_df, abs_wi_df, rel_wi_df
