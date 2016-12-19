import unittest

from alfpy import ncd
from alfpy.utils import distmatrix

from . import utils


class Test(unittest.TestCase, utils.ModulesCommonTest):

    def __init__(self, *args, **kwargs):
        super(Test, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()

    def test_complexity1(self):
        seq = 'AACGTACCATTGAACGTACCGTAGG'
        c = ncd.complexity(seq)
        self.assertEqual(c, 26)

    def test_complexity2(self):
        seq = 'MFTDNAKIIIVQLNASVEINCTRPNNNTR'
        c = ncd.complexity(seq)
        self.assertEqual(c, 37)

    def test_complexities(self):
        dist = ncd.Distance(self.pep_records)
        exp = [
            ((0,), 63.0), ((0, 1), 77.0), ((0, 2), 85.0),
            ((0, 3), 70.0), ((1,), 60.0), ((1, 2), 78.0),
            ((1, 3), 65.0), ((2,), 61.0), ((2, 3), 66.0),
            ((3,), 37.0)
        ]
        self.assertEqual(exp, sorted(dist._complexity.items()))

    def test_distance(self):
        dist = ncd.Distance(self.pep_records)
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.2698413 0.3809524 0.5238095",
            "seq2       0.2698413 0.0000000 0.2950820 0.4666667",
            "seq3       0.3809524 0.2950820 0.0000000 0.4754098",
            "seq4       0.5238095 0.4666667 0.4754098 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))


if __name__ == '__main__':
    unittest.main()
