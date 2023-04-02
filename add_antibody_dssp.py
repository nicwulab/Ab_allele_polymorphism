from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import pandas as pd
from Bio.PDB.DSSP import DSSP
import sys
from itertools import combinations
from itertools import permutations




def main(pdb_dir, mkdssp_dir):
    file_format = '.pdb'
    p = PDBParser()
    err = []
    result_df = pd.read_csv('results/epitope_identification_with_antibody_only_ddgs.csv', index_col=0)
    d = {}
    for index, row_data in result_df.iterrows():
        try:
            pdb_id = row_data['pdb_id'].lower()

            row_chain_id = row_data['chain_id']

            row_dssp_location = row_data['dssp_location']
            row_dssp_location_digit = ''
            row_dssp_location_alpha = ''

            for s in row_dssp_location:
                if s.isdigit():
                    row_dssp_location_digit += s
                else:
                    row_dssp_location_alpha += s
            row_dssp = str(row_dssp_location_digit)+row_dssp_location_alpha.strip()
            dssp_location_digit = int(row_dssp_location_digit)
            structure = p.get_structure(pdb_id, pdb_dir + pdb_id + file_format)
            model = structure[0]

            dssp = DSSP(model, pdb_dir + pdb_id + file_format, mkdssp_dir)

            for i in range(len(list(dssp.keys()))):

                chain_id = str(dssp.keys()[i][0])
                res_subres_id = str(dssp.keys()[i][1][1]) + str(dssp.keys()[i][1][2])
                d[(pdb_id, chain_id, res_subres_id)] = dssp[list(dssp.keys())[i]][3]

                if chain_id == row_chain_id and res_subres_id.strip() == row_dssp:
                    result_df.loc[index, 'dssp_asa'] = dssp[list(dssp.keys())[i]][3]

        except:
            err.append(index)
    print(err)
    result_df.to_csv('results/add_dssp_area.csv')
    import json
    my_dict = {str(key): value for key, value in d.items()}

    with open("dssp_asa.txt", "w") as f:
        json.dump(my_dict, f)

if __name__ == "__main__":
    pdb_dir = sys.argv[1]
    mkdssp_dir = sys.argv[2]
    main(pdb_dir, mkdssp_dir)


