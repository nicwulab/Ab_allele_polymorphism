import unittest

import part_2 as part_2

class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, False)  # add assertion here

    def test_is_equal(self):
        anarci = [('H1', 'Q'), ('H2', 'V'), ('H3', 'Q'), ('H4', 'L'), ('H5', 'V'), ('H6', 'Q'), ('H7', 'S'),
                  ('H8', 'G'), ('H9', 'A'), ('H11', 'E'), ('H12', 'V'), ('H13', 'K'), ('H14', 'K'), ('H15', 'P'),
                  ('H16', 'G'), ('H17', 'S')]

        dssp = [('H1', 'X'), ('H1', 'X'), ('H2', 'V'), ('H3', 'Q'), ('H1', 'X'), ('H1', 'X'), ('H4', 'L'), ('H5', 'V'),
                ('H6', 'Q'), ('H7', 'S'), ('H8', 'G'), ('H9', 'X'), ('H10', 'E'), ('H11', 'V'), ('H12', 'K'),
                ('H13', 'K'), ('H14', 'P'), ('H15', 'G'), ('H16', 'S')]

        self.assertEqual(part_2.equal(anarci, dssp, 0, 0, []), False)



if __name__ == '__main__':
    unittest.main()
