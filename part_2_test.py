import unittest
from part_2 import equal, get_anarci_location_for_part_1_result, igv_csv_to_df, pdb_seq_csv_to_df, dict_to_geneid, \
    get_igv_seq_df_dssp_chain, get_positions_interact_with_antigen_from_part_1_result, trim_dssp_chain_and_part_1_result
import pandas as pd
import utils.helper as helper

class MyTestCase(unittest.TestCase):

    def test_is_equal(self):
        anarci = [('H17', 'S')]
        dssp = [('H16', 'S')]

        self.assertEqual(equal(anarci, dssp, 0, 0, []), True)

    def test_is_equal_true(self):
        anarci = [('H1', 'Q'), ('H2', 'V'), ('H3', 'Q'), ('H4', 'L'), ('H5', 'V'), ('H6', 'Q'), ('H7', 'S'),
                  ('H8', 'G'), ('H9', 'A'), ('H11', 'E'), ('H12', 'V'), ('H13', 'K'), ('H14', 'K'), ('H15', 'P'),
                  ('H16', 'G'), ('H17', 'S')]
        dssp = [('H1', 'X'), ('H1', 'X'), ('H2', 'V'), ('H3', 'Q'), ('H1', 'X'), ('H1', 'X'), ('H4', 'L'), ('H5', 'V'),
                ('H6', 'Q'), ('H7', 'S'), ('H8', 'G'), ('H9', 'X'), ('H10', 'E'), ('H11', 'V'), ('H12', 'K'),
                ('H13', 'K'), ('H14', 'P'), ('H15', 'G'), ('H16', 'S')]

        self.assertEqual(equal(anarci, dssp, 0, 0, []), True)

    def test_is_equal_false(self):
        anarci = [('H1', 'Q'), ('H2', 'V'), ('H3', 'Q'), ('H4', 'L'), ('H5', 'V'), ('H6', 'Q'), ('H7', 'S'),
                  ('H8', 'G'), ('H9', 'A'), ('H11', 'E'), ('H12', 'V'), ('H13', 'K'), ('H14', 'K'), ('H15', 'P'),
                  ('H16', 'G'), ('H17', 'S')]
        dssp = [('H1', 'A'), ('H1', 'X'), ('H2', 'V'), ('H3', 'Q'), ('H1', 'X'), ('H1', 'X'), ('H4', 'L'), ('H5', 'V'),
                ('H6', 'Q'), ('H7', 'S'), ('H8', 'G'), ('H9', 'X'), ('H10', 'E'), ('H11', 'V'), ('H12', 'K'),
                ('H13', 'K'), ('H14', 'P'), ('H15', 'G'), ('H16', 'S')]

        self.assertEqual(equal(anarci, dssp, 0, 0, []), False)

    def test_is_equal_false_case_2(self):
        anarci = [('H1', 'Q'), ('H2', 'V')]
        dssp = [('H1', 'A'), ('H1', 'X')]

        self.assertEqual(equal(anarci, dssp, 0, 0, []), False)

    def test_is_equal_true_case_2(self):
        anarci = [('H1', 'Q'), ('H2', 'V'), ('H3', 'Q'), ('H4', 'L'), ('H5', 'V'), ('H6', 'Q'), ('H7', 'S'),
                  ('H8', 'G'), ('H9', 'A'), ('H11', 'E'), ('H12', 'V'), ('H13', 'K'), ('H14', 'K'), ('H15', 'P'),
                  ('H16', 'G'), ('H17', 'S')]
        dssp = [('H1', 'X'), ('H1', 'X'), ('H2', 'X'), ('H2', 'X'), ('H3', 'Q'), ('H1', 'X'), ('H1', 'X'), ('H4', 'L'),
                ('H5', 'V'), ('H6', 'Q'), ('H7', 'S'), ('H8', 'G'), ('H9', 'X'), ('H10', 'E'), ('H11', 'V'), ('H12', 'K'),
                ('H13', 'K'), ('H14', 'P'), ('H15', 'G'), ('H16', 'S')]
        self.assertEqual(equal(anarci, dssp, 0, 0, []), True)

    def test_get_anarci_location_dict(self):
        pdb_chain = '6PHD:H'
        dssp_chain_locations = [('H1', 'A'), ('H2', 'V'), ('H3', 'Q'), ('H5', 'V'), ('H16', 'S')]
        anarci_row_data = pd.DataFrame({'Id': '6PHD:H', '2': 'A', '3': 'V', '4': 'Q', '5': 'V', '6': 'S'}, index=[0])

        expected_dict = {'H1': 'H2', 'H2': 'H3', 'H3': 'H4', 'H16': 'H6'}
        self.assertEqual(get_anarci_location_for_part_1_result(pdb_chain, dssp_chain_locations, anarci_row_data),
                         expected_dict)

    def test_get_anarci_location_dict_anarci_shorter_than_dssp(self):
        pdb_chain = '6PHD:H'
        dssp_chain_locations = \
            [('H1', 'A'), ('H2', 'V'), ('H3', 'Q'), ('H5', 'V'), ('H16', 'S'), ('H17', 'S'), ('H18', 'S'), ('H19', 'S')]
        anarci_row_data = pd.DataFrame({'Id': '6PHD:H', '2': 'A', '3': 'V', '4': 'Q', '5': 'V', '6': 'S'}, index=[0])

        expected_dict = {'H1': 'H2', 'H2': 'H3', 'H3': 'H4', 'H16': 'H6'}
        self.assertEqual(get_anarci_location_for_part_1_result(pdb_chain, dssp_chain_locations, anarci_row_data), expected_dict)

    def test_get_anarci_location_dict_anarci_contains_empty_values(self):
        pdb_chain = '6PHD:H'
        dssp_chain_locations = \
            [('H1', 'A'), ('H3', 'Q'), ('H5', 'V'), ('H16', 'S'), ('H17', 'S'), ('H18', 'S'), ('H19', 'S')]
        anarci_row_data = pd.DataFrame({'Id': '6PHD:H', '2': 'A', '3': '-', '4': 'Q', '5': 'V', '6': 'S'}, index=[0])

        expected_dict = {'H1': 'H2', 'H3': 'H4', 'H16': 'H6'}
        self.assertEqual(get_anarci_location_for_part_1_result(pdb_chain, dssp_chain_locations, anarci_row_data),
                         expected_dict)

    def test_get_anarci_location_dict_dssp_contains_x(self):
        pdb_chain = '6PHD:H'
        dssp_chain_locations = [('H1', 'A'), ('H2', 'X'), ('H3', 'Q'), ('H5', 'V'), ('H16', 'S')]
        anarci_row_data = pd.DataFrame({'Id': '6PHD:H', '2': 'A', '3': 'V', '4': 'Q', '5': 'V', '6': 'S'}, index=[0])

        expected_dict = {'H1': 'H2', 'H2': 'H3', 'H3': 'H4', 'H16': 'H6'}
        self.assertEqual(get_anarci_location_for_part_1_result(pdb_chain, dssp_chain_locations, anarci_row_data),
                         expected_dict)

    def test_dict_to_geneid(self):
        list_of_dict = [{'Sequence ID': '3SE8:L',
                         'Raw Sequence': 'EIVLTQSPGILSLSP', 'Sequence Length': 208,
                         'Domain Classification': 'imgt', 'Hits': [{'gene': 'IGKV3-20*01 unnamed protein product',
                                                                    'bit_score': 131.0, 'e_value': 3e-42},
                                                                   {'gene': 'IGKV3D-20*02 unnamed protein product',
                                                                    'bit_score': 130.0, 'e_value': 7e-42},
                                                                   {'gene': 'IGKV3D-11*01 unnamed protein product',
                                                                    'bit_score': 128.0, 'e_value': 8e-41}]},
                        {'Sequence ID': '3MXW:H', 'Raw Sequence': 'QVQLQQSGPELVRP',
                                                        'Sequence Length': 220, 'Domain Classification': 'imgt',
                                                        'Hits': [{'gene': 'IGHV1-2*02 unnamed protein product',
                                                                  'bit_score': 129.0, 'e_value': 4e-41},
                                                                 {'gene': 'IGHV1-2*06 unnamed protein product',
                                                                  'bit_score': 129.0, 'e_value': 4e-41},
                                                                 {'gene': 'IGHV1-2*05 unnamed protein product',
                                                                  'bit_score': 127.0, 'e_value': 2e-40}]}]
        expected_list = [['3SE8:L', 'EIVLTQSPGILSLSP', 'IGKV3-20'],
                         ['3MXW:H', 'QVQLQQSGPELVRP', 'IGHV1-2']]

        self.assertEqual(dict_to_geneid(list_of_dict), expected_list)

    def test_get_igv_seq_df_dssp_chain(self):
        gene_id = 'IGHV1-2*02'
        dssp_pdb_location = [[('H1 ', 'E'), ('H2 ', 'V'), ('H3 ', 'Q'), ('H4 ', 'L'), ('H5 ', 'V'), ('H6 ', 'E')],
                             [('L3 ', 'E'), ('L4 ', 'L'), ('L5 ', 'T'), ('L6 ', 'Q'), ('L7 ', 'E'), ('L8 ', 'T'),
                              ('L9 ', 'G'), ('L11 ', 'V'), ('L12 ', 'S'), ('L13 ', 'V'), ('L14 ', 'A'), ('L15 ', 'L'),
                              ('L16 ', 'G'), ('L17 ', 'Q'), ('L18 ', 'T'), ('L19 ', 'V'), ('L20 ', 'T'), ('L21 ', 'I'),
                              ('L22 ', 'T'), ('L23 ', 'C'), ('L24 ', 'Q'), ('L25 ', 'G'), ('L26 ', 'D'), ('L27 ', 'S'),
                              ('L28 ', 'L'), ('L29 ', 'R'), ('L30 ', 'S'), ('L31 ', 'H'), ('L32 ', 'Y'), ('L33 ', 'A')]]
        igv_H_df = igv_csv_to_df('data/anarci_igv_output.csv_H.csv')
        igv_KL_df = igv_csv_to_df('data/anarci_igv_output.csv_KL.csv')

        pdb_seq_H_df = pdb_seq_csv_to_df('data/anarci_pdb_output.csv_H.csv')
        pdb_seq_KL_df = pdb_seq_csv_to_df('data/anarci_pdb_output.csv_KL.csv')

        _, _, chain_locations = get_igv_seq_df_dssp_chain(gene_id, dssp_pdb_location, igv_H_df, igv_KL_df, pdb_seq_H_df, pdb_seq_KL_df)
        self.assertEqual(chain_locations, [('H1', 'E'), ('H2', 'V'), ('H3', 'Q'), ('H4', 'L'), ('H5', 'V'), ('H6', 'E')])

    def test_get_positions_interact_with_antigen_from_part_1_result(self):
        part_1_result_dict = {'2b1a': [('L30 ', 'N'), ('L31 ', 'N'), ('L32 ', 'Y'), ('L50 ', 'R'), ('L51 ', 'N'),
                                       ('L66 ', 'K'), ('L91 ', 'W'), ('L92 ', 'D'), ('L93 ', 'D'), ('L95A', 'G'),
                                       ('L95B', 'G'), ('L95C', 'P'), ('H31 ', 'D'), ('H32 ', 'Y'), ('H33 ', 'W'),
                                       ('H50 ', 'I'), ('H52 ', 'Y'), ('H54 ', 'D'), ('H56 ', 'D'), ('H57 ', 'S'),
                                       ('H58 ', 'R'), ('H95 ', 'L'), ('H98 ', 'D'), ('H99 ', 'Y'), ('H100 ', 'E'),
                                       ('H100A', 'D'), ('H100B', 'S'), ('H100C', 'G'), ('H100D', 'A'), ('H100E', 'D')]}
        pdb = '2B1A'
        chain_id = 'L'
        expected_list = [('L30', 'N'), ('L31', 'N'), ('L32', 'Y'), ('L50', 'R'), ('L51', 'N'), ('L66', 'K'), ('L91', 'W'),
                         ('L92', 'D'), ('L93', 'D'), ('L95A', 'G'), ('L95B', 'G'), ('L95C', 'P')]
        self.assertEqual(get_positions_interact_with_antigen_from_part_1_result(part_1_result_dict, pdb.lower(), chain_id),
                         expected_list)

    def test_trim_dssp_chain_and_part_1_result_trim_end(self):
        dssp_pdb_location = [('H1', 'E'), ('H2', 'V'), ('H3', 'Q'), ('H4', 'L'), ('H5', 'V'), ('H6', 'E')]
        start_index = 0
        end_index = 2
        position_to_index = {'H1': 0, 'H2': 1, 'H3': 2, 'H4': 3, 'H5': 4, 'H6': 5}
        positions_interact_with_antigen_tuples = [('H1', 'E'), ('H2', 'V'), ('H3', 'Q')]

        adjusted_chain_locations, positions_interact_with_antigen_tuples = \
            trim_dssp_chain_and_part_1_result(dssp_pdb_location, start_index, end_index, position_to_index,
                                              positions_interact_with_antigen_tuples)

        self.assertEqual(adjusted_chain_locations, [('H1', 'E'), ('H2', 'V'), ('H3', 'Q')])
        self.assertEqual(positions_interact_with_antigen_tuples, [('H1', 'E'), ('H2', 'V'), ('H3', 'Q')])

    def test_trim_dssp_chain_and_part_1_result_trim_start(self):
        dssp_pdb_location = [('H1', 'E'), ('H2', 'V'), ('H3', 'Q'), ('H4', 'L'), ('H5', 'V'), ('H6', 'E')]
        start_index = 2
        end_index = 5
        position_to_index = {'H1': 0, 'H2': 1, 'H3': 2, 'H4': 3, 'H5': 4, 'H6': 5}
        positions_interact_with_antigen_tuples = [('H1', 'E'), ('H2', 'V'), ('H3', 'Q'), ('H4', 'L'), ('H5', 'V'),
                                                  ('H6', 'E')]

        adjusted_chain_locations, positions_interact_with_antigen_tuples = \
            trim_dssp_chain_and_part_1_result(dssp_pdb_location, start_index, end_index, position_to_index,
                                              positions_interact_with_antigen_tuples)

        self.assertEqual(adjusted_chain_locations, [('H3', 'Q'), ('H4', 'L'), ('H5', 'V'), ('H6', 'E')])
        self.assertEqual(positions_interact_with_antigen_tuples, [('H3', 'Q'), ('H4', 'L'), ('H5', 'V'), ('H6', 'E')])

    def test_trim_dssp_chain_and_part_1_result_trim_start_and_end(self):
        dssp_pdb_location = [('H1', 'E'), ('H2', 'V'), ('H3', 'Q'), ('H4', 'L'), ('H5', 'V'), ('H6', 'E')]
        start_index = 2
        end_index = 4
        position_to_index = {'H1': 0, 'H2': 1, 'H3': 2, 'H4': 3, 'H5': 4, 'H6': 5}
        positions_interact_with_antigen_tuples = [('H1', 'E'), ('H2', 'V'), ('H3', 'Q'), ('H4', 'L'), ('H5', 'V'),
                                                  ('H6', 'E')]

        adjusted_chain_locations, positions_interact_with_antigen_tuples = \
            trim_dssp_chain_and_part_1_result(dssp_pdb_location, start_index, end_index, position_to_index,
                                              positions_interact_with_antigen_tuples)

        self.assertEqual(adjusted_chain_locations, [('H3', 'Q'), ('H4', 'L'), ('H5', 'V')])
        self.assertEqual(positions_interact_with_antigen_tuples, [('H3', 'Q'), ('H4', 'L'), ('H5', 'V')])

    def test_trim_dssp_chain_and_part_1_result_trim_no_trim(self):
        dssp_pdb_location = [('H1', 'E'), ('H2', 'V'), ('H3', 'Q'), ('H4', 'L')]
        start_index = 0
        end_index = 3
        position_to_index = {'H1': 0, 'H2': 1, 'H3': 2, 'H4': 3}
        positions_interact_with_antigen_tuples = [('H1', 'E'), ('H2', 'V')]

        adjusted_chain_locations, positions_interact_with_antigen_tuples = \
            trim_dssp_chain_and_part_1_result(dssp_pdb_location, start_index, end_index, position_to_index,
                                              positions_interact_with_antigen_tuples)

        self.assertEqual(adjusted_chain_locations, [('H1', 'E'), ('H2', 'V'), ('H3', 'Q'), ('H4', 'L')])
        self.assertEqual(positions_interact_with_antigen_tuples, [('H1', 'E'), ('H2', 'V')])


if __name__ == '__main__':

    unittest.main()
