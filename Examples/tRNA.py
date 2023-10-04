# First we are going to filter the records and save them in another file.
# By this, we will also see the filter_records function for advanced users

from CodonU.analyzer import filter_reference  # for filtering records
from Bio.SeqIO import parse, write  # for reading the existing .fasta file and writing filtered file

handle = 'Nucleotide/Staphylococcus_aureus.fasta'
records = parse(handle, 'fasta')
filtered_records = filter_reference(records, min_len_threshold=200, _type='nuc')

new_file_path = 'Nucleotide/Staphylococcus_aureus_tRNA_200.fasta'  # tRNA defines the use, and 200 is the minimum length
with open(new_file_path, 'w') as out_file:
    write(filtered_records, out_file, 'fasta')
    out_file.close()

# Now we will first retrieve the tRNA gene copy number from GtRNAdb

from CodonU.analyzer import get_anticodon_count_dict

anti_dict = get_anticodon_count_dict(
    url='http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/bacteria/Stap_aure_aureus_Newman/',
    database='GtRNAdb'
)
# just for seeing what's inside the dict, execute the following block
# for codon, num in anti_dict.items():
#     print(f'{codon}: {num}')

# Lastly to calculate gtAI, we will import the function and give necessary arguments

from CodonU.analyzer import calculate_gtai

tai_df, abs_wi_df, rel_wi_df = calculate_gtai(
    handle='Nucleotide/Staphylococcus_aureus_tRNA_200.fasta',
    anticodon_dict=anti_dict,
    genetic_code_num=11,
    generation_num=20
)

# for seeing the contents, execute below block of codes
# print("The tAI values are in a Dataframe, and is as below:")
# print(tai_df)
# print('############################################################')
# print("The absolute values for weights are in a Dataframe, and is as below:")
# print(abs_wi_df)
# print('############################################################')
# print("The relative values for weights are in a Dataframe, and is as below:")
# print(rel_wi_df)
# print('############################################################')
