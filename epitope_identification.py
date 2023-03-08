from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import pandas as pd
from Bio.PDB.DSSP import DSSP
import sys
from itertools import combinations


def get_list_antigen_loc(pdb_dir, mkdssp_dir):
    result_df = pd.read_csv('results/compute_ddG/compute_ddG_result.csv', index_col=[0])

    file_format = '.pdb'
    p = PDBParser()
    errors = []
    emptys = []
    res_dic = {}

    for index, row_data in result_df.iterrows():

        if row_data['pdb_id'].lower() in res_dic:
            result_df.loc[index, 'antigen_loc_list'] = res_dic[row_data['pdb_id'].lower()]
            continue
        try:
            pdb_id = row_data['pdb_id'].lower()
            chain_id = row_data['chain_id']
            antigen_chain_id = row_data['antigen_chain']

            dssp_location = row_data['dssp_location']
            dssp_location_digit = ''
            dssp_location_alpha = ''

            for s in dssp_location:
                if s.isdigit():
                    dssp_location_digit += s
                else:
                    dssp_location_alpha += s

            dssp_location_digit = int(dssp_location_digit)
            structure = p.get_structure(pdb_id, pdb_dir + pdb_id + file_format)
            model = structure[0]

            dssp = DSSP(model, pdb_dir + pdb_id + file_format, mkdssp_dir)

            antigen_chain_id_to_asa = pd.DataFrame(columns=['chain_res_subres_id', 'asa', 'amino_acid'])
            # print(dssp.keys())
            for i in range(len(list(dssp.keys()))):
                chain_res_subres_id = str(dssp.keys()[i][0]) + str(dssp.keys()[i][1][1]) + str(dssp.keys()[i][1][2])
                if str(dssp.keys()[i][0]) in antigen_chain_id:
                    new_row = {'chain_res_subres_id': chain_res_subres_id, 'asa': dssp[list(dssp.keys())[i]][3],
                               'amino_acid': dssp[list(dssp.keys())[i]][1]}
                    antigen_chain_id_to_asa = antigen_chain_id_to_asa.append(new_row, ignore_index=True)

            # Remove all chains except light chain and heavy chain
            input_file = pdb_dir + pdb_id + file_format
            output_file = "without_target.pdb"

            structure_vl_vh = p.get_structure("pdb", input_file)

            for m in structure_vl_vh:
                for chain in list(m):
                    if chain.id not in antigen_chain_id:
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

            df_remove_res_id_to_asa = pd.DataFrame(columns=['chain_res_subres_id', 'asa_without_target', 'amino_acid'])

            for i in range(len(list(dssp_removed.keys()))):
                chain_res_subres_id = str(dssp_removed.keys()[i][0]) + str(dssp_removed.keys()[i][1][1]) + str(dssp_removed.keys()[i][1][2])
                if str(dssp_removed.keys()[i][0]) in antigen_chain_id:
                    new_row = {'chain_res_subres_id': chain_res_subres_id,
                               'asa_without_target': dssp_removed[list(dssp_removed.keys())[i]][3],
                               'amino_acid': dssp_removed[list(dssp_removed.keys())[i]][1]}
                    df_remove_res_id_to_asa = df_remove_res_id_to_asa.append(new_row, ignore_index=True)

            df2 = pd.merge(antigen_chain_id_to_asa, df_remove_res_id_to_asa, on=['chain_res_subres_id'])
            df3 = df2[df2['asa'] < df2['asa_without_target']]
            antigen_chain_res_subres_id_list = df3['chain_res_subres_id']

            antigen_chain_res_subres_id_list = [s.strip() for s in antigen_chain_res_subres_id_list.tolist()]

            str_antigen_chain_res_subres_id_list = ','.join(antigen_chain_res_subres_id_list)
            res_dic[row_data['pdb_id'].lower()] = str_antigen_chain_res_subres_id_list
            result_df.loc[index, 'antigen_loc_list'] = str_antigen_chain_res_subres_id_list
        except:
            errors.append((index, row_data['pdb_id'].lower()))

    result_df.to_csv('epitope_identification.csv')
    print(errors)
    print(emptys)

def set_map(set1, set2):
    import itertools

    intersection_lenghts = [[0 for j in range(len(set1))] for i in range(len(set1))]

    for i in range(len(set1)):
        for j in range(len(set2)):
            intersection_lenghts[i][j] = len(set1[i].intersection(set2[j]))
    a = [i for i in range(len(set1))]

    permutations = itertools.permutations(a)
    mappings = []
    for p in permutations:
        mapping = tuple(zip(a, p))
        mappings.append(mapping)

    max_intersection_length = 0
    max_intersection_combo = ()

    for t in mappings:
        local_intersection_length = 0
        for mapping in t:
            local_intersection_length += intersection_lenghts[mapping[0]][mapping[1]]
        if local_intersection_length > max_intersection_length:
            max_intersection_length = local_intersection_length
            max_intersection_combo = t

    set_length = 0
    for s in set1:
        set_length += len(s)
    for s in set2:
        set_length += len(s)

    matching_percent = max_intersection_length / (set_length / 2)

    return matching_percent


def get_similarity_for_diff_chain_name(group_id, group_df, result_df, diff_chain):
    first_chain_name_list = group_df['antigen_chain'].unique()[0].split('|')
    first_chain_name_list = [c.strip() for c in first_chain_name_list]

    antigen_loc_lists = {}
    for index, row_data in group_df.iterrows():
        antigen_loc_list_df = row_data['antigen_loc_list'].split(',')
        antigen_chain_list = row_data['antigen_chain'].split('|')
        antigen_chain_list = [c.strip() for c in antigen_chain_list]

        antigen_loc_lists[index] = {}
        for i in range(len(antigen_chain_list)):
            chain_loc = []
            for antigen_loc in antigen_loc_list_df:
                if antigen_chain_list[i] in antigen_loc:
                    chain_loc.append(antigen_loc[1:])

            antigen_loc_lists[index][i] = chain_loc

    group_ids = {}
    next_group_id = 1

    for list1, list2 in combinations(antigen_loc_lists.keys(), 2):
        positions1 = antigen_loc_lists[list1]
        positions1 = positions1.values()
        set1 = []

        for p in positions1:
            set1.append(set(p))

        positions2 = antigen_loc_lists[list2]
        positions2 = positions2.values()
        set2 = []
        for p in positions2:
            set2.append(set(p))

        matching_percent = set_map(set1, set2)

        if matching_percent >= 0.8:
            if list1 in group_ids and list2 in group_ids:
                # Merge the two groups into a single group
                group_id = group_ids[list1]
                for key, value in group_ids.items():
                    if value == group_ids[list2]:
                        group_ids[key] = group_id
            elif list1 in group_ids:
                # Add list2 to the same group as list1
                group_id = group_ids[list1]
                group_ids[list2] = group_id
            elif list2 in group_ids:
                # Add list1 to the same group as list2
                group_id = group_ids[list2]
                group_ids[list1] = group_id
            else:
                # Create a new group
                group_ids[list1] = next_group_id
                group_ids[list2] = next_group_id
                next_group_id += 1

    for l in antigen_loc_lists.keys():
        if l not in group_ids:
            group_ids[l] = next_group_id
            next_group_id += 1

    for index, val in group_ids.items():
        result_df.loc[index, 'groupID'] = str(result_df.loc[index, 'groupID']) + '-' + str(val)
        group_df.loc[index, 'groupID'] = str(group_df.loc[index, 'groupID']) + '-' + str(val)


def get_group_id():
    species_id_df = pd.read_csv('species_ID.csv')
    species_id_df = species_id_df.drop('Unnamed: 2', axis=1)
    species_id_df = species_id_df.drop('Unnamed: 3', axis=1)
    species_id_df = species_id_df.drop('Unnamed: 4', axis=1)

    species_to_id_dict = dict(zip(species_id_df.species, species_id_df.ID))

    result_df = pd.read_csv('epitope_identification_sorted.csv', index_col=[0])

    for index, row_data in result_df.iterrows():
        species_list = row_data['antigen_species']
        if not isinstance(species_list, str):
            result_df.loc[index, 'antigen_species_id'] = None
            continue
        species_list = species_list.split('|')
        species_list = [s.strip() for s in species_list]
        id_list = []
        for s in species_list:
            id = species_to_id_dict[s]
            id_list.append(id)

        result_df.loc[index, 'antigen_species_id'] = ','.join(id_list)

    def assign_group_id(df):
        return df.groupby(['antigen_species_id', 'gene']).ngroup()

    result_df['groupID'] = assign_group_id(result_df)

    result_df.to_csv('test.csv')
    grouped_df = result_df.groupby('groupID')

    more_than_one_chain_group_ids = []
    for group_id, group_df in grouped_df:
        if group_id == -1:
            continue

        pdb_list = group_df['pdb_id'].unique()
        if len(pdb_list) == 1:
            continue

        diff_chain = []
        if ',' in group_df['antigen_species_id'].unique()[0]:
            diff_chain.append(group_id)
            print("here")
            get_similarity_for_diff_chain_name(group_id, group_df, result_df, diff_chain)
            continue

        antigen_loc_lists = {}

        for index, row_data in group_df.iterrows():
            antigen_loc_list = row_data['antigen_loc_list'].split(',')
            antigen_loc_list = [s[1:] for s in antigen_loc_list]
            antigen_loc_lists[index] = antigen_loc_list

        group_ids = {}
        next_group_id = 1
        threshold = 0.8

        for list1, list2 in combinations(antigen_loc_lists.keys(), 2):

            set1 = set(antigen_loc_lists[list1])
            set2 = set(antigen_loc_lists[list2])

            intersection = set1.intersection(set2)
            num_elements = (len(set1) + len(set2)) / 2

            similarity = len(intersection) / num_elements

            if similarity >= threshold:
                if list1 in group_ids and list2 in group_ids:
                    # Merge the two groups into a single group
                    group_id = group_ids[list1]
                    for key, value in group_ids.items():
                        if value == group_ids[list2]:
                            group_ids[key] = group_id
                elif list1 in group_ids:
                    # Add list2 to the same group as list1
                    group_id = group_ids[list1]
                    group_ids[list2] = group_id

                elif list2 in group_ids:
                    # Add list1 to the same group as list2
                    group_id = group_ids[list2]
                    group_ids[list1] = group_id
                else:
                    # Create a new group
                    group_ids[list1] = next_group_id
                    group_ids[list2] = next_group_id
                    next_group_id += 1

        for l in antigen_loc_lists.keys():
            if l not in group_ids:
                group_ids[l] = next_group_id
                next_group_id += 1

        for index, val in group_ids.items():
            result_df.loc[index, 'groupID'] = str(result_df.loc[index, 'groupID']) + '-' + str(val)
            group_df.loc[index, 'groupID'] = str(group_df.loc[index, 'groupID']) + '-' + str(val)

    unique_groupID = result_df['groupID'].unique().tolist()

    groupID_dict = {}

    for i in range(len(unique_groupID)):
        if unique_groupID[i] == -1:
            continue
        groupID_dict[unique_groupID[i]] = str(int(i + 1))

    # print(groupID_dict)
    result_df['groupID'] = result_df['groupID'].map(groupID_dict)

    result_df.to_csv('test4.csv')

def main(pdb_dir, mkdssp_dir):
    # get_list_antigen_loc(pdb_dir, mkdssp_dir)

    get_group_id()


if __name__ == "__main__":
    pdb_dir = sys.argv[1]
    mkdssp_dir = sys.argv[2]
    main(pdb_dir, mkdssp_dir)