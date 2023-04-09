from os.path import join, abspath
from CodonU.analyzer import calculate_cai, calculate_enc, calculate_rscu, calculate_cbi, calculate_gravy, \
    calculate_aromaticity
from CodonU.file_handler import make_dir
from CodonU.file_handler.internal_comp import is_file, is_file_empty
from CodonU.cua_errors import UnsupportedType
from Bio.SeqUtils import seq3


def _print_msg(_type: str, output, tot_gene: int):
    """
    Prints the upper message
    """
    print('The report is generated with the help of package CodonU created by Souradipto C.',
          file=output)
    print('_' * 80, file=output)
    print(f'Total gene count: {tot_gene}', file=output)
    print('Value Ranges:', file=output)
    if _type == 'nuc':
        print('    RSCU:', file=output)
        print('        - > 1: There is a positive codon usage bias', file=output)
        print('        - < 1: There is a negative codon usage bias', file=output)
        print('    CAI:', file=output)
        print('        - 0: The usage is highly biased', file=output)
        print('        - 1: The usage is completely random', file=output)
        print('        - nan: There are no synonymous codons', file=output)
        print('    CBI:', file=output)
        print('        - 0: The usage is highly biased', file=output)
        print('        - 1: The usage is completely random', file=output)
        print('        - nan: There are no synonymous codons', file=output)
        print('    ENc:', file=output)
        print('        - 61: The usage is completely random', file=output)
        print('        - 20: The usage is highly biased', file=output)
    elif _type == 'aa':
        print('    GRAVY:', file=output)
        print('        - > 0: The protein is relatively hydrophobic', file=output)
        print('        - < 0: The protein is relatively hydrophilic', file=output)
        print('        - nan: There are no synonymous codons', file=output)
        print('    Aromaticity:', file=output)
        print('        - 0: The protein consists 0 aromatic residues', file=output)
        print('        - > 0: The protein consists aromatic residues', file=output)
        print('        - Higher the value is, higher are the number of aromatics residues in the protein', file=output)


def generate_report(handle: str, _type: str, genetic_code_num: int, min_len_threshold: int,
                    res_folder_path: str = 'Report'):
    """
    Generate the report for given sequence **[best for gene analysis]**

    For nucleotide sequence, this generates reports of:
        - RSCU
        - CAI
        - CBI
        - ENc

    For protein sequence, this generates reports of:
        - GRAVY score
        - Aromaticity score


    **NOTE** Possible types are
        - nuc: For nucleotide sequence
        - aa: For protein sequence

    :param handle: Handle to the file, or the filename as a string
    :param _type: Type of the sequence [nuc or aa]
    :param genetic_code_num: Genetic table number for codon table
    :param min_len_threshold: Minimum length of sequence to be considered as gene
    :param res_folder_path: The path of folder where the file will be saved
    :return:
    """
    supported_types = ['nuc', 'aa']
    if _type not in supported_types:
        raise UnsupportedType
    make_dir(res_folder_path)
    identifier = handle.split('/')[-1].split('.')[0]
    report_file_name = f'report_{identifier}_nuc.txt' if _type == 'nuc' else f'report_{identifier}_aa.txt'
    report_file_path = join(res_folder_path, report_file_name)
    gene_analysis = True
    if not is_file(report_file_path) or is_file_empty(report_file_path):
        if _type == 'nuc':
            print('Calculating RSCU, please be patient, this may take some time.')
            rscu_dict = calculate_rscu(handle, genetic_code_num, min_len_threshold, gene_analysis)
            print('Calculating CAI, please be patient, this may take some time.')
            cai_dict = calculate_cai(handle, genetic_code_num, min_len_threshold, gene_analysis)
            tot_gene = len(cai_dict.keys())
            print('Calculating CBI, please be patient, this may take some time.')
            cbi_dict = calculate_cbi(handle, genetic_code_num, min_len_threshold, gene_analysis)
            print('Calculating ENc, please be patient, this may take some time.')
            enc_dict = calculate_enc(handle, genetic_code_num, min_len_threshold, gene_analysis)
            with open(report_file_path, 'w') as output:
                _print_msg(_type, output, tot_gene)
                print('\nRSCU Details', file=output)
                for key, _dict in rscu_dict.items():
                    print(f'    {key}:', file=output)
                    for codon, val in _dict.items():
                        print(f'        {codon}: {val}', file=output)
                print('\nCAI Details', file=output)
                for key, _dict in cai_dict.items():
                    print(f'    {key}:', file=output)
                    for codon, val in _dict.items():
                        print(f'        {codon}: {val}', file=output)
                print('\nCBI Details', file=output)
                for key, _dict in cbi_dict.items():
                    print(f'    {key}:', file=output)
                    for aa, val in _dict.items():
                        print(f'        {seq3(aa)}: {val[0]} (Most used codon {val[1]})', file=output)
                print('\nENc Details', file=output)
                for key, val in enc_dict.items():
                    print(f'    {key}: {val}', file=output)
        elif _type == 'aa':
            print('Calculating GRAVY score, please be patient, this may take some time.')
            gravy_dict = calculate_gravy(handle, min_len_threshold, gene_analysis)
            tot_gene = len(gravy_dict.keys())
            print('Calculating aromaticity, please be patient, this may take some time.')
            aroma_dict = calculate_aromaticity(handle, min_len_threshold, gene_analysis)
            with open(report_file_path, 'w') as output:
                _print_msg(_type, output, tot_gene)
                print('\nGRAVY Details', file=output)
                for key, val in gravy_dict.items():
                    print(f'    {key}: {val}', file=output)
                print('\nAromaticity Details', file=output)
                for key, val in aroma_dict.items():
                    print(f'    {key}: {val}', file=output)
        print(f'\nThe report can be found at {abspath(report_file_path)}')
