import pandas as pd
import utils.helper as helper

not_match_count = 0
not_match_list = []
processed_list = []
no_chain_location = []
empty_border_list = []


def dict_to_geneid(list_of_dict):
    list_of_id_rawseq_geneid = []
    for d in list_of_dict:
        id = d['Sequence ID']
        seq = d['Raw Sequence']
        gene_id_dict = d['Hits'][0]
        gene = gene_id_dict['gene']
        gene_id = gene.split('*')[0]
        list_of_id_rawseq_geneid.append([id, seq, gene_id])

    return list_of_id_rawseq_geneid


def igv_csv_to_df(filename):
    df = pd.read_csv(filename)
    df[['Id', 'allele']] = df['Id'].str.split('*', expand=True)
    df = df.drop('domain_no', axis=1)
    df = df.drop('hmm_species', axis=1)
    df = df.drop('e-value', axis=1)
    df = df.drop('identity_species', axis=1)
    df = df.drop('v_gene', axis=1)
    df = df.drop('v_identity', axis=1)
    df = df.drop('j_gene', axis=1)
    df = df.drop('j_identity', axis=1)
    df = df.drop('seqstart_index', axis=1)
    df = df.drop('score', axis=1)
    df = df.drop('seqend_index', axis=1)
    df = df.drop('chain_type', axis=1)
    return df


def pdb_seq_csv_to_df(filename):
    df = pd.read_csv(filename)
    df = df.drop('domain_no', axis=1)
    df = df.drop('hmm_species', axis=1)
    df = df.drop('e-value', axis=1)
    df = df.drop('identity_species', axis=1)
    df = df.drop('v_gene', axis=1)
    df = df.drop('v_identity', axis=1)
    df = df.drop('j_gene', axis=1)
    df = df.drop('j_identity', axis=1)
    df = df.drop('score', axis=1)
    df = df.drop('chain_type', axis=1)
    return df


def get_anarci_location_for_part_1_result(pdb_chain, chain_locations, row_data):
    pdb, chainid = pdb_chain.split(':')

    row_data = row_data.iloc[0]
    row_tuples = [(col, val) for col, val in row_data.items()]
    row_tuples_without_id = row_tuples[1:]
    row_tuples_with_chain_id = [(chainid + s[0], s[1]) for s in row_tuples_without_id]

    anarci_row_tuples_without_empty_slots = [s for s in row_tuples_with_chain_id if s[1] != '-']

    dssp_location_to_anarci = dict()

    result = []
    is_two_equal = equal(anarci_row_tuples_without_empty_slots, chain_locations, 0, 0, result)
    if not is_two_equal:
        global not_match_count
        global not_match_list
        not_match_count += 1
        not_match_list.append(pdb_chain)
        return {'error': 'not equal'}

    for i in result:
        del chain_locations[i]

    for i in range(len(anarci_row_tuples_without_empty_slots)):
        anarci_pos, anarci_name = anarci_row_tuples_without_empty_slots[i]
        dssp_pos, dssp_name = chain_locations[i]

        if anarci_pos != dssp_pos:
            dssp_location_to_anarci[dssp_pos] = anarci_pos

    return dssp_location_to_anarci


def equal(anarci, dssp, i, j, result):
    if i == len(anarci) or j == len(dssp):
        return True

    if dssp[j][1] == anarci[i][1]:
        return equal(anarci, dssp, i + 1, j + 1, result)

    if dssp[j][1] == 'X':
        if equal(anarci, dssp, i + 1, j + 1, result):
            return True
        else:
            try_stay = equal(anarci, dssp, i, j + 1, result)
            if try_stay:
                result.append(j)
                return True
    return False


def get_igv_seq_df_dssp_chain(gene_id, dssp_pdb_location, igv_H_df, igv_KL_df, pdb_seq_H_df, pdb_seq_KL_df):
    hchain_locations, lchain_locations = dssp_pdb_location
    igv = None
    pdb_seq_df = None
    chain_locations = None

    if gene_id[:3] == 'IGH':
        igv = igv_H_df
        pdb_seq_df = pdb_seq_H_df
        chain_locations = hchain_locations
    elif gene_id[:3] == 'IGK' or gene_id[:3] == 'IGL':
        igv = igv_KL_df
        pdb_seq_df = pdb_seq_KL_df
        chain_locations = lchain_locations

    chain_locations = [(s[0].strip(), s[1]) for s in chain_locations]

    return igv, pdb_seq_df, chain_locations


def get_positions_interact_with_antigen_from_part_1_result(part_1_result_dict, pdb, chain_id):
    positions_interact_with_antigen_raw = part_1_result_dict[pdb]
    positions_interact_with_antigen_tuples = [(t[0].strip(), t[1]) for t in
                                              positions_interact_with_antigen_raw]
    positions_interact_with_antigen_tuples = [s for s in positions_interact_with_antigen_tuples if
                                              s[0][0] == chain_id]

    return positions_interact_with_antigen_tuples


def trim_dssp_chain_and_part_1_result(chain_locations, start_index, end_index, position_to_index,
                                      positions_interact_with_antigen_tuples):

    adjusted_chain_locations = chain_locations[start_index:end_index + 1]
    discarded_locations = chain_locations[:start_index]
    discarded_locations_end = chain_locations[end_index + 1:]

    # if end_index < len(chain_locations), need to remove tuples in the positions_interact_with_antigen_tuples
    # list that are after end_index
    if len(discarded_locations_end) > 0:
        new_positions_interact_with_antigen_tuples = []
        last_discarded_position_str = discarded_locations_end[0][0]
        for t in positions_interact_with_antigen_tuples:
            position, name = t
            if position_to_index[position] < position_to_index[last_discarded_position_str]:
                new_positions_interact_with_antigen_tuples.append(t)
        positions_interact_with_antigen_tuples = new_positions_interact_with_antigen_tuples

    # if start_index > 0, need to remove tuples in the positions_interact_with_antigen_tuples list that are before
    # start_index
    if start_index != 0:
        new_positions_interact_with_antigen_tuples = []
        last_discarded_position_str = discarded_locations[-1][0]
        for t in positions_interact_with_antigen_tuples:
            position, name = t
            if position_to_index[position] > position_to_index[last_discarded_position_str]:
                new_positions_interact_with_antigen_tuples.append(t)
        positions_interact_with_antigen_tuples = new_positions_interact_with_antigen_tuples

    return adjusted_chain_locations, positions_interact_with_antigen_tuples


def compare_each_id_to_gene_id(part_1_result_dict, list_of_id_rawseq_geneid, igv_H_df, igv_KL_df,
                               dssp_pdb_location_dict, pdb_seq_H_df, pdb_seq_KL_df, pdb_chain_compound_name_df):
    errors = []
    result_df = pd.DataFrame(
        columns=['pdb_id', 'chain_id', 'location', 'gene', 'amino_acid_original', 'list_amino_acid_variants', 'compound'])

    for pdb_chain, _, gene_id in list_of_id_rawseq_geneid:
        pdb, chain_id = pdb_chain.split(':')
        if pdb.lower() not in part_1_result_dict:
            errors.append(pdb.lower())
            continue

        if len(part_1_result_dict[pdb.lower()]) == 0:
            global empty_border_list
            empty_border_list.append(pdb_chain)
            continue

        dssp_pdb_location = dssp_pdb_location_dict[pdb.lower()]
        igv, pdb_seq_df, chain_locations = get_igv_seq_df_dssp_chain(gene_id, dssp_pdb_location, igv_H_df, igv_KL_df,
                                                                     pdb_seq_H_df, pdb_seq_KL_df)

        pdb_seq_row = pdb_seq_df[pdb_seq_df['Id'] == pdb_chain]

        for index, row_data in pdb_seq_row.iterrows():

            positions_interact_with_antigen_tuples = \
                get_positions_interact_with_antigen_from_part_1_result(part_1_result_dict, pdb.lower(), chain_id)

            row_data = row_data.to_frame()  # series to frame
            row_data = row_data.T
            start_index = int(row_data['seqstart_index'])
            end_index = int(row_data['seqend_index'])

            row_data = row_data.drop('seqstart_index', axis=1)
            row_data = row_data.drop('seqend_index', axis=1)

            # build map: location_name to index in the list
            position_to_index = dict()
            for i in range(len(chain_locations)):
                position_to_index[chain_locations[i][0]] = i

            adjusted_chain_locations, positions_interact_with_antigen_tuples = \
                trim_dssp_chain_and_part_1_result(chain_locations, start_index, end_index, position_to_index, positions_interact_with_antigen_tuples)

            # after trimming chain_locations, map dssp location to anarci location
            dssp_location_to_anarci = get_anarci_location_for_part_1_result(pdb_chain, adjusted_chain_locations,
                                                                            row_data)
            if 'error' in dssp_location_to_anarci:
                continue

            # create map for interacting amino acid locations to ANARCI position
            for i in range(len(positions_interact_with_antigen_tuples)):
                dssp_position, amino_acid_name = positions_interact_with_antigen_tuples[i]
                if dssp_position in dssp_location_to_anarci:
                    positions_interact_with_antigen_tuples[i] = (dssp_location_to_anarci[dssp_position], amino_acid_name)

            # remove chain id from the list of positions_interact_with_antigen_tuples (i.e. H12 -> 12)
            positions_interact_with_antigen_dict = dict(positions_interact_with_antigen_tuples)
            positions_interact_with_antigen = positions_interact_with_antigen_dict.keys()
            positions_without_chainid = [s[1:] for s in positions_interact_with_antigen]

            # for each column of the igv dataframe where igv['Id'] == gene_id, compare the the rows, if there are more
            # than 2 unique values then check if they are in the list of positions_interact_with_antigen,
            # if yes, write to result dataframe
            gene_id_df = igv[igv['Id'] == gene_id]
            compound_name = pdb_chain_compound_name_df[pdb_chain_compound_name_df['pdb'] == pdb.lower()]['compound']
            for col in gene_id_df.columns:
                if len(gene_id_df[col].unique()) > 1:
                    if len(gene_id_df[col].unique()) == 2 and '-' in gene_id_df[col].unique():
                        continue
                    if col == 'allele':
                        continue
                    if col in positions_without_chainid:
                        list_variants = gene_id_df[col].unique()
                        if '-' in list_variants:
                            list_variants = list_variants[list_variants != '-']
                        new_row = {'pdb_id': pdb,
                                   'chain_id': chain_id,
                                   'location': col,
                                   'gene': gene_id,
                                   'amino_acid_original': positions_interact_with_antigen_dict[chain_id + col],
                                   'list_amino_acid_variants': list_variants,
                                   'compound': compound_name}
                        result_df = result_df.append(new_row, ignore_index=True)

            global processed_list
            processed_list.append(pdb_chain)

    print("errors", errors)
    print("errors count", len(errors))
    result_df.to_csv('part_2_result_2_test.csv')
    return


def main():
    list_of_dict = helper.parse_pyir_output('results/part1/pdb_sequences.json')
    test_list = ['7WEE:H','7JWB:D']
    # test_list = ['6IEA:L', '6M3B:C', '4HPO:L', '6IUT:L', '7SD5:L', '7SJP:H', '4YE4:L', '6VJN:L', '3H0T:A', '5HHV:L',
    #              '5HHX:L', '8DWA:L', '7M7B:L', '4XC3:L', '4XCF:H', '4RIS:L', '3JCB:B', '3JCC:B', '2H32:A', '5D70:L',
    #              '7TTY:L', '7TTX:L', '7TTM:L', '7JVA:L', '7D03:H', '7D6I:C', '6W7S:L', '7WO7:B', '7WOG:B', '7RP2:H',
    #              '7N4M:L', '5SY8:L']
    list_of_dict = [s for s in list_of_dict]
    # list_of_dict = [s for s in list_of_dict if s['Sequence ID'] in test_list]

    print(len(list_of_dict))
    list_of_id_rawseq_geneid = dict_to_geneid(list_of_dict)

    pdb_chain_compound_name_df = pd.read_table('utils/pdb_with_only_one_chain_pairing_with_compound_name.tsv')

    igv_H_df = igv_csv_to_df('data/anarci_igv_output.csv_H.csv')
    igv_KL_df = igv_csv_to_df('data/anarci_igv_output.csv_KL.csv')

    pdb_seq_H_df = pdb_seq_csv_to_df('data/anarci_pdb_output.csv_H.csv')
    pdb_seq_KL_df = pdb_seq_csv_to_df('data/anarci_pdb_output.csv_KL.csv')

    part_1_result_dict = helper.parse_result_with_tuple('results/part1/part_1_result_with_amino_acid_name.txt')
    dssp_pdb_location_dict = helper.parse_result_with_tuple('results/part1/dssp_pdb_seq_location.txt')

    compare_each_id_to_gene_id(part_1_result_dict, list_of_id_rawseq_geneid, igv_H_df, igv_KL_df,
                               dssp_pdb_location_dict, pdb_seq_H_df, pdb_seq_KL_df, pdb_chain_compound_name_df)

    print('not match list', not_match_list)
    print('not match count', not_match_count)
    print('processed list', processed_list)
    print('processed count', len(processed_list))
    print('no chain', no_chain_location)
    print('no chain count', len(no_chain_location))
    print('empty_border_list', empty_border_list)
    print('empty_border_list count', len(empty_border_list))


if __name__ == "__main__":
    main()
