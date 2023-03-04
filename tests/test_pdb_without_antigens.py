from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import pandas as pd
from Bio.PDB.DSSP import DSSP
from Bio import SeqIO
from Bio.SeqUtils import seq1
import json
import ast

def parse_summary_file(filename):
    df = pd.read_table(filename)
    df = df.drop_duplicates(subset=['pdb'], keep=False)
    df_HL = df[['pdb', 'Hchain', 'Lchain', 'antigen_chain']]
    print(df_HL)
    df3 = df_HL[df_HL['antigen_chain'].isnull()]
    return df3

def parse_result_file(filename):
    file1 = open(filename, 'r')
    string = file1.readline()
    json_acceptable_string = string.replace("'", "\"")
    d = json.loads(json_acceptable_string)
    return d

def parse_result_with_tuple(filename):
    # Open the file for reading
    with open(filename, 'r') as f:
        # Read the contents of the file as a string
        contents = f.read()

    # Use ast.literal_eval to convert the string into a dictionary
    return ast.literal_eval(contents)

unique_pdbs_no_antigen_df = parse_summary_file('../data/20230217_0084705_summary.tsv')
print(unique_pdbs_no_antigen_df)

result_dict = (parse_result_file('../results/pdb_to_paratope/pdb_to_paratope_result.txt'))
print(len(result_dict.keys()))
result_dict_with_name = (parse_result_with_tuple(
    '../results/pdb_to_paratope/pdb_to_paratope_result_with_amino_acid_name.txt'))
print(len(result_dict_with_name.keys()))
count_elements = 0
for k, v in result_dict_with_name.items():
    count_elements += len(v)
print(count_elements)

unique_pdbs_no_antigen_list = unique_pdbs_no_antigen_df['pdb'].tolist()
print(unique_pdbs_no_antigen_list)
print(len(unique_pdbs_no_antigen_list))

count = 0
for u in unique_pdbs_no_antigen_list:
    if u not in result_dict:
        continue
    if result_dict[u] != []:
        print(u,':', result_dict[u])
        count += 1
print(count)
#


