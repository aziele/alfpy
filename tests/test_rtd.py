import unittest

from alfpy import rtd
from alfpy import word_pattern
from alfpy.utils import distmatrix

from . import utils


class Test(unittest.TestCase, utils.ModulesCommonTest):
    """Shared methods and tests for creating RTD vectors."""

    def __init__(self, *args, **kwargs):
        super(Test, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()
        self.pep_2mer_pos = word_pattern.create(
            self.pep_records.seq_list, 2, True)

    def test_calc_rtd(self):
        seq = 'CTACACAACTTTGCGGGTAGCCGGAAACATTGTGAATGCGGTGAACA'
        apos = [i for i, nt in enumerate(seq) if nt == 'A']
        val = rtd.calc_rtd(apos, 1)
        exp = (3.3846153846153846, 3.1510306381944679)
        self.assertEqual(val, exp)

    def test_create_vector(self):
        vector = rtd.create_vector(self.pep_records.count, self.pep_2mer_pos)
        exp = (self.pep_records.count, len(self.pep_2mer_pos.pat_list)*2)
        self.assertEqual(vector.shape, exp)

    def test_distance(self):
        vector = rtd.create_vector(self.pep_records.count, self.pep_2mer_pos)
        dist = rtd.Distance(vector, 'google')
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.4892241 0.6034483 0.9310345",
            "seq2       0.4892241 0.0000000 0.3673469 0.8802817",
            "seq3       0.6034483 0.3673469 0.0000000 0.8843537",
            "seq4       0.9310345 0.8802817 0.8843537 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

if __name__ == '__main__':
    unittest.main()
