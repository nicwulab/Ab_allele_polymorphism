from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import pandas as pd
from Bio.PDB.DSSP import DSSP
from Bio import SeqIO


def get_dssp_position_amino_acid_dict(pdb, hchain, lchain):
    p = PDBParser()
    dir = '/Users/natalieso/Downloads/20230217_0084705/'
    file_format = '.pdb'
    structure = p.get_structure(pdb, dir + pdb + file_format)

    model = structure[0]
    dssp = DSSP(model, dir + pdb + file_format, '/Users/natalieso/Downloads/dssp-3.1.4/mkdssp')

    hchain_list = []
    lchain_list = []
    for i in range(len(dssp.keys())):
        k = dssp.keys()[i]
        v = dssp[k]
        chain = k[0]
        if chain == lchain:
            lchain_list.append((k[0] + str(k[1][1]) + k[1][2], v[1]))
        if chain == hchain:
            hchain_list.append((k[0] + str(k[1][1]) + k[1][2], v[1]))

    return [hchain_list, lchain_list]

def parse_summary_file(filename):
    df = pd.read_table(filename)
    df = df.drop_duplicates(subset=['pdb'], keep=False)
    df_HL = df[['pdb', 'Hchain', 'Lchain']]
    df_HL.to_csv('pdb_with_only_one_chain_pairing.tsv', sep="\t")

    return df_HL

unique_pdbs_df = parse_summary_file('../data/20230217_0084705_summary.tsv')
pdb_list = unique_pdbs_df['pdb'].tolist()

unique_pdbs_df.set_index("pdb", inplace=True)

dir = '/Users/natalieso/Downloads/20230217_0084705/'

result_dict = dict()
errors = list()

for pdb in pdb_list:
    pdb_row = unique_pdbs_df.loc[pdb]
    lchain = pdb_row['Lchain']
    hchain = pdb_row['Hchain']

    try:
        lists = get_dssp_position_amino_acid_dict(pdb, hchain, lchain)
        result_dict[pdb] = lists
    except:
        errors.append(pdb)

print(len(result_dict.keys()))
print(result_dict)
with open('dssp_pdb_seq_location.txt', 'w') as data:
    data.write(str(result_dict))

print(errors)