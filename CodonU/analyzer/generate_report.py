from os.path import join
from Bio.SeqIO import parse
from CodonU.analyzer import calculate_cai, calculate_enc, calculate_rscu, calculate_cbi, calculate_gravy, \
    calculate_aromaticity
from CodonU.file_handler import make_dir
from CodonU.file_handler.internal_comp import is_file, is_file_empty


# TODO implement detailed report writing function for protein sequence
# TODO details about _type

def generate_report(handle: str, _type: str, genetic_code_num: int, min_len_threshold: int,
                    res_folder_path: str = 'Report', gene_analysis: bool = False):
    make_dir(res_folder_path)
    identifier = handle.split('/')[-1].split('.')[0]
    report_file_name = f'report_{identifier}.txt'
    report_file_path = join(res_folder_path, report_file_name)
    if not is_file(report_file_path) or is_file_empty(report_file_path):
        if _type == 'gene':
            print('Calculating CAI')
            records = parse(handle, 'fasta')
            cai_dict = calculate_cai(records, genetic_code_num, min_len_threshold, gene_analysis)
            print('Calculating RSCU')
            records = parse(handle, 'fasta')
            rscu_dict = calculate_rscu(records, genetic_code_num, min_len_threshold, gene_analysis)
            print('Calculating ENc')
            records = parse(handle, 'fasta')
            enc_dict = calculate_enc(records, genetic_code_num, min_len_threshold, gene_analysis)
            with open(report_file_path, 'w') as output:
                print('The report is generated with the help of package CodonU created by Souradipto C.',
                      file=output)
                if gene_analysis:
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
                    print('\nENc Details', file=output)
                    for key, val in enc_dict.items():
                        print(f'    {key}: {val}', file=output)
                else:
                    print('\nRSCU Details', file=output)
                    for key, val in rscu_dict.items():
                        print(f'    {key}:{val}', file=output)
                    print('\nCAI Details', file=output)
                    for key, val in cai_dict.items():
                        print(f'    {key}:{val}', file=output)
                    print('\nENc Details', file=output)
                    print(f'ENc value: {enc_dict}', file=output)
        elif _type == 'aa':
            pass


if __name__ == '__main__':
    handle = '../../Results/Nucleotide/temp.fasta'
    generate_report(handle, 'gene', 11, 200, gene_analysis=False)
