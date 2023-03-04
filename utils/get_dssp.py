from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import pandas as pd
from Bio.PDB.DSSP import DSSP
from Bio import SeqIO
import sys

def get_dssp_position_amino_acid_dict(pdb, hchain, lchain, pdb_dir, mkdssp_dir):
    p = PDBParser()
    dir = pdb_dir
    file_format = '.pdb'
    structure = p.get_structure(pdb, dir + pdb + file_format)

    model = structure[0]
    dssp = DSSP(model, dir + pdb + file_format, mkdssp_dir)

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


def main(pdb_dir, summary_file, mkdssp_dir):
    # get positions of interacting amino acids for each pdb with only one chain pairings
    unique_pdbs_df = parse_summary_file(summary_file)
    pdb_list = unique_pdbs_df['pdb'].tolist()

    unique_pdbs_df.set_index("pdb", inplace=True)

    dir = pdb_dir

    result_dict = dict()
    errors = list()

    for pdb in pdb_list:
        pdb_row = unique_pdbs_df.loc[pdb]
        lchain = pdb_row['Lchain']
        hchain = pdb_row['Hchain']

        try:
            lists = get_dssp_position_amino_acid_dict(pdb, hchain, lchain, pdb_dir, mkdssp_dir)
            result_dict[pdb] = lists
        except:
            errors.append(pdb)

    with open('results/pdb_to_paratope/dssp_pdb_seq_location.txt', 'w') as data:
        data.write(str(result_dict))

    print(errors)

if __name__ == "__main__":
    pdb_dir = sys.argv[1]
    summary_file = sys.argv[2]
    mkdssp_dir = sys.argv[3]
    main(pdb_dir, summary_file, mkdssp_dir)