from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import pandas as pd
from Bio.PDB.DSSP import DSSP
import sys


def parse_summary_file(filename):
    df = pd.read_table(filename)
    df = df.drop_duplicates(subset=['pdb'], keep=False)
    df_HL = df[['pdb', 'Hchain', 'Lchain', 'compound']]
    df_HL.to_csv('utils/pdb_with_only_one_chain_pairing.tsv', sep="\t")

    return df_HL

def parse_pdb_file(pdb_name, df, dir, mkdssp_dir):
    file_format = '.pdb'

    p = PDBParser()
    structure = p.get_structure(pdb_name, dir+pdb_name+file_format)

    pdb_row = df.loc[pdb_name]
    lchain = pdb_row['Lchain']
    hchain = pdb_row['Hchain']

    model = structure[0]

    dssp = DSSP(model, dir+pdb_name+file_format, mkdssp_dir)

    df_all_chain_id_to_asa = pd.DataFrame(columns=['chain_res_subres_id', 'asa', 'amino_acid'])

    for i in range(len(list(dssp.keys()))):

        chain_res_subres_id = str(dssp.keys()[i][0]) + str(dssp.keys()[i][1][1]) + str(dssp.keys()[i][1][2])
        new_row = {'chain_res_subres_id':chain_res_subres_id, 'asa': dssp[list(dssp.keys())[i]][3], 'amino_acid': dssp[list(dssp.keys())[i]][1]}
        df_all_chain_id_to_asa = df_all_chain_id_to_asa.append(new_row, ignore_index=True)

    # Remove all chains except light chain and heavy chain
    input_file = dir+pdb_name+file_format
    output_file = "VL_VH.pdb"

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

    model_removed = structure_vl_vh[0]

    dssp_removed = DSSP(model_removed, output_file, mkdssp_dir)

    df_remove_chain_id_to_asa = pd.DataFrame(columns=['chain_res_subres_id', 'asa_without_antigen', 'amino_acid'])

    for i in range(len(list(dssp_removed.keys()))):
        chain_res_subres_id = str(dssp_removed.keys()[i][0]) + str(dssp_removed.keys()[i][1][1]) + str(dssp_removed.keys()[i][1][2])
        new_row = {'chain_res_subres_id': chain_res_subres_id, 'asa_without_antigen': dssp_removed[list(dssp_removed.keys())[i]][3], 'amino_acid': dssp_removed[list(dssp_removed.keys())[i]][1]}
        df_remove_chain_id_to_asa = df_remove_chain_id_to_asa.append(new_row, ignore_index=True)

    df2 = pd.merge(df_all_chain_id_to_asa, df_remove_chain_id_to_asa, on=['chain_res_subres_id'])
    if not df2['amino_acid_x'].equals(df2['amino_acid_y']):
        print("amino acid name not equal")
        return

    df3 = df2[df2['asa'] < df2['asa_without_antigen']]

    return df3

def get_positions(pdb_dir, summary_file, mkdssp_dir):
    # parse summary file to filter out duplicates
    # returns a df of pdbs with Hchain and Lchain
    unique_pdbs_df = parse_summary_file(summary_file)
    pdb_list = unique_pdbs_df['pdb'].tolist()

    unique_pdbs_df.set_index("pdb", inplace=True)

    dir = pdb_dir

    result_dict = dict()

    errors = list()
    for pdb in pdb_list:
        try:
            result_df = parse_pdb_file(pdb, unique_pdbs_df, dir, mkdssp_dir)
            result_df['pos_name'] = list(zip(result_df.chain_res_subres_id, result_df.amino_acid_y))

            result_dict[pdb] = result_df['pos_name'].tolist()
        except:
            errors.append(pdb)

    with open('pdb_to_paratope_result_with_amino_acid_name.txt', 'w') as data:
        data.write(str(result_dict))

    print(errors)

def main(pdb_dir, summary_file, mkdssp_dir):
    # get positions of interacting amino acids for each pdb with only one chain pairings
    get_positions(pdb_dir, summary_file, mkdssp_dir)

if __name__ == "__main__":
    pdb_dir = sys.argv[1]
    summary_file = sys.argv[2]
    mkdssp_dir = sys.argv[3]
    main(pdb_dir, summary_file, mkdssp_dir)

