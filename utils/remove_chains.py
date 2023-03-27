from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import pandas as pd
from Bio.PDB.DSSP import DSSP
from Bio import SeqIO
import sys


def parse_summary_file(filename):
    df = pd.read_table(filename)
    df = df.drop_duplicates(subset=['pdb'], keep=False)
    df_HL = df[['pdb', 'Hchain', 'Lchain']]
    df_HL.to_csv('pdb_with_only_one_chain_pairing.tsv', sep="\t")

    return df_HL


def main(pdb_dir, summary_file, output_dir):
    # get positions of interacting amino acids for each pdb with only one chain pairings
    unique_pdbs_df = parse_summary_file(summary_file)
    pdb_list = unique_pdbs_df['pdb'].tolist()

    unique_pdbs_df.set_index("pdb", inplace=True)

    dir = pdb_dir

    for pdb_name in pdb_list:
        file_format = '.pdb'

        p = PDBParser()

        pdb_row = unique_pdbs_df.loc[pdb_name]
        lchain = pdb_row['Lchain']
        hchain = pdb_row['Hchain']


        # Remove all chains except light chain and heavy chain
        input_file = dir + pdb_name + file_format
        output_file = output_dir+"/"+pdb_name+".pdb"

        structure_vl_vh = p.get_structure("pdb", input_file)

        for m in structure_vl_vh:
            for chain in list(m):
                if chain.id == lchain or chain.id == hchain:
                    continue
                m.detach_child(chain.id)

        # Create a PDB writer object and write out the modified structure
        writer = PDBIO()
        writer.set_structure(structure_vl_vh)
        writer.save(output_file)
        f = open(output_file)
        text = f.read()
        f.close()
        f = open(output_file, 'w')
        f.write("CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n")
        f.write(text)
        f.close()


if __name__ == "__main__":
    pdb_dir = sys.argv[1]
    summary_file = sys.argv[2]
    output_dir = sys.argv[3]
    main(pdb_dir, summary_file, output_dir)