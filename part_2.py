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
    row_tuples_with_chain_id = [(chainid+s[0], s[1]) for s in row_tuples_without_id]

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


def compare_each_id_to_gene_id(part_1_result_dict, list_of_id_rawseq_geneid, igv_H_df, igv_KL_df, dssp_pdb_location_dict):
    errors = []
    result_df = pd.DataFrame(columns=['pdb_id', 'chain_id', 'location', 'gene', 'amino_acid_original', 'list_amino_acid_variants'])

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
        else:
            print("more igv?")

        chain_locations = [(s[0].strip(), s[1]) for s in chain_locations]

        pdb_seq_row = pdb_seq_df[pdb_seq_df['Id'] == pdb_chain]

        for index, row_data in pdb_seq_row.iterrows():

            positions_interact_with_antigen_raw = part_1_result_dict[pdb.lower()]
            positions_interact_with_antigen_tuples = [(t[0].strip(), t[1]) for t in positions_interact_with_antigen_raw]
            positions_interact_with_antigen_tuples = [s for s in positions_interact_with_antigen_tuples if
                                                      s[0][0] == chain_id]

            row_data = row_data.to_frame()
            row_data = row_data.T
            start_index = int(row_data['seqstart_index'])
            end_index = int(row_data['seqend_index'])

            row_data = row_data.drop('seqstart_index', axis=1)
            row_data = row_data.drop('seqend_index', axis=1)

            adjusted_chain_locations = chain_locations[start_index:end_index+1]
            discarded_locations = chain_locations[:start_index]
            discarded_locations_end = chain_locations[end_index+1:]

            position_to_index = dict()
            for i in range(len(chain_locations)):
                position_to_index[chain_locations[i][0]] = i
            if len(discarded_locations_end) > 0:
                new_positions_interact_with_antigen_tuples = []
                last_discarded_position_str = discarded_locations_end[0][0]
                for t in positions_interact_with_antigen_tuples:
                    position, name = t
                    if position_to_index[position] < position_to_index[last_discarded_position_str]:
                        new_positions_interact_with_antigen_tuples.append(t)
                positions_interact_with_antigen_tuples = new_positions_interact_with_antigen_tuples

            if start_index != 0:
                new_positions_interact_with_antigen_tuples = []
                last_discarded_position_str = discarded_locations[-1][0]

                for t in positions_interact_with_antigen_tuples:
                    position, name = t
                    if position_to_index[position] > position_to_index[last_discarded_position_str]:
                        new_positions_interact_with_antigen_tuples.append(t)

                positions_interact_with_antigen_tuples = new_positions_interact_with_antigen_tuples

            dssp_location_to_anarci = get_anarci_location_for_part_1_result(pdb_chain, adjusted_chain_locations, row_data)
            if 'error' in dssp_location_to_anarci:
                continue

            # create map for interacting amino acid locations to ANARCI position
            for i in range(len(positions_interact_with_antigen_tuples)):
                dssp_position, amino_acid_name = positions_interact_with_antigen_tuples[i]
                if dssp_position in dssp_location_to_anarci:
                    positions_interact_with_antigen_tuples[i] = (dssp_location_to_anarci[dssp_position], amino_acid_name)

            positions_interact_with_antigen_dict = dict(positions_interact_with_antigen_tuples)
            positions_interact_with_antigen = positions_interact_with_antigen_dict.keys()
            positions_without_chainid = [s[1:] for s in positions_interact_with_antigen]
            gene_id_df = igv[igv['Id'] == gene_id]

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
                                   'amino_acid_original': positions_interact_with_antigen_dict[chain_id+col],
                                   'list_amino_acid_variants': list_variants}
                        result_df = result_df.append(new_row, ignore_index=True)

            global processed_list
            processed_list.append(pdb_chain)

    print("errors", errors)
    print("errors count", len(errors))
    result_df.to_csv('part_2_result_2_test.csv')
    return


list_of_dict = helper.parse_pyir_output('results/part1/pdb_sequences.json')
# print(list_of_dict)
# test_list = ['7WEE:H', '6OZ4:H', '7KRA:I', '7BYR:H', '6WO3:H', '6IEA:L', '7NIW:C', '7CYV:H', '7BEM:L', '4XVJ:H', '4XVJ:L', '7RAN:E', '6M3B:C', '7U9P:H', '7U9O:H', '6NOU:A', '4HPY:L', '7O4Y:L', '7AJ6:H', '4HPO:L', '6CW2:B', '4JPW:H', '1AR2:A', '6IUT:L', '7SD5:L', '6MLK:L', '7D5U:D', '7SJP:H', '4YE4:L', '4LRN:L', '6VJN:L', '6G8R:A', '3H0T:A', '6OZ2:H', '5HHV:H', '5HHV:L', '5HHX:L', '3KYK:H', '6CBV:H', '6CBP:A', '6E65:L', '8DWA:L', '5ERW:A', '7M7B:L', '7KBM:A', '7N3C:L', '7KBO:A', '6VBO:L', '4XC3:L', '6U1N:L', '4XCF:H', '4R26:L', '4RIS:L', '6QB3:X', '3QOT:L', '6DSI:A', '3JCB:B', '3JCC:B', '6TCR:H', '2H32:A', '4OB5:L', '4FQ1:L', '4FQ2:L', '4FQC:L', '7TI6:L', '6QFC:B', '6URH:H', '5I8O:H', '5WT9:L', '7L77:L', '5VKD:H', '6QB4:X', '7LO8:L', '7LO7:H', '6VUP:A', '6VUG:C', '6VUG:D', '7JXC:L', '5D70:L', '4UT7:L', '7RDJ:A', '7ARN:H', '7ARN:L', '7RDM:A', '7RDK:A', '6BKB:H', '7BZ5:H', '7BZ5:L', '7EPU:A', '5Y11:A', '4HWB:H', '7TTY:L', '7TTX:L', '7TTM:L', '5E08:L', '7JVA:L', '3KR3:L', '7F1G:D', '4HQQ:L', '7O52:L', '5IR3:A', '1X9Q:A', '5VAG:B', '6TCQ:H', '4JY5:L', '4JY4:A', '7D03:H', '3BDY:H', '7DX4:H', '3DVI:A', '7XP4:N', '7XP6:N', '6U8D:L', '7D6I:C', '6PZG:L', '7U9W:H', '6YXF:H', '7D6Y:B', '2B4C:H', '8DCE:H', '1VGE:H', '7DAA:L', '7KTX:I', '3BJ9:1', '5DRW:B', '6UKJ:H', '6W7S:L', '5I18:L', '6EHW:B', '6EHV:X', '6EHX:B', '7WO7:B', '7WOG:B', '7WOQ:D', '7RP2:H', '7A3T:L', '7N4M:L', '6J9O:H', '5SY8:L', '7CWO:L', '6YX9:L', '6YXG:H', '7ZCF:E']
# test_list= ['6YXG:H', '4KQ3:H', '7TI6:L', '6QFC:B', '7WEE:H', '7JWB:D']
# test_list = ['7WEE:H','7JWB:D']
# test_list = ['6IEA:L', '6M3B:C', '4HPO:L', '6IUT:L', '7SD5:L', '7SJP:H', '4YE4:L', '6VJN:L', '3H0T:A', '5HHV:L', '5HHX:L', '8DWA:L', '7M7B:L', '4XC3:L', '4XCF:H', '4RIS:L', '3JCB:B', '3JCC:B', '2H32:A', '5D70:L', '7TTY:L', '7TTX:L', '7TTM:L', '7JVA:L', '7D03:H', '7D6I:C', '6W7S:L', '7WO7:B', '7WOG:B', '7RP2:H', '7N4M:L', '5SY8:L']
# test_list = ['8DCE:H']
list_of_dict = [s for s in list_of_dict]
# list_of_dict = [s for s in list_of_dict if s['Sequence ID'] in test_list]
print(len(list_of_dict))
list_of_id_rawseq_geneid = dict_to_geneid(list_of_dict)

igv_H_df = igv_csv_to_df('data/anarci_igv_output.csv_H.csv')
igv_KL_df = igv_csv_to_df('data/anarci_igv_output.csv_KL.csv')

pdb_seq_H_df = pdb_seq_csv_to_df('data/anarci_pdb_output.csv_H.csv')
pdb_seq_KL_df = pdb_seq_csv_to_df('data/anarci_pdb_output.csv_KL.csv')

part_1_result_dict = helper.parse_result_with_tuple('results/part1/part_1_result_with_amino_acid_name.txt')
dssp_pdb_location_dict = helper.parse_result_with_tuple('results/part1/dssp_pdb_seq_location.txt')

compare_each_id_to_gene_id(part_1_result_dict, list_of_id_rawseq_geneid, igv_H_df, igv_KL_df, dssp_pdb_location_dict)

print('not match list', not_match_list)
print('not match count', not_match_count)
print('processed list', processed_list)
print('processed count', len(processed_list))
print('no chain', no_chain_location)
print('no chain count', len(no_chain_location))
print('empty_border_list', empty_border_list)
print('empty_border_list count', len(empty_border_list))




