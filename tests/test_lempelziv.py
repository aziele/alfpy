import unittest

from alfpy import lempelziv
from alfpy.utils import distmatrix

from . import utils


class VectorTest(unittest.TestCase, utils.ModulesCommonTest):

    def __init__(self, *args, **kwargs):
        super(VectorTest, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()

    def test_complexity(self):
        seq = 'MFTDNAKIIIVQLNASVEINCTRPNNNTR'
        c = lempelziv.complexity(seq)
        self.assertEqual(c, 19)

    def test_complexity1(self):
        seq = 'MFTDNAKIIIVQLNASVEINCTRPNNNTR'
        c = lempelziv.complexity1(seq)
        self.assertEqual(c, 20)

    def test_complexities(self):
        dist = lempelziv.Distance(self.pep_records)
        exp = [((0,), 40), ((0, 1), 47), ((0, 2), 53),
               ((0, 3), 43), ((1,), 38), ((1, 0), 47),
               ((1, 2), 47), ((1, 3), 41), ((2,), 35),
               ((2, 0), 50), ((2, 1), 45), ((2, 3), 37),
               ((3,), 19), ((3, 0), 39), ((3, 1), 37),
               ((3, 2), 36)]
        self.assertEqual(sorted(dist._complexity.items()), exp)


class DistanceTest(unittest.TestCase, utils.ModulesCommonTest):

    def __init__(self, *args, **kwargs):
        super(DistanceTest, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()
        self.dist = lempelziv.Distance(self.pep_records, 'd')

    def test_distance_d(self):
        matrix = distmatrix.create(self.pep_records.id_list, self.dist)
        exp = [
            "   4",
            "seq1       0 9 15 20",
            "seq2       9 0 10 18",
            "seq3       15 10 0 17",
            "seq4       20 18 17 0"
        ]
        self.assertEqual(matrix.format(decimal_places=0), "\n".join(exp))

    def test_distance_d_star(self):
        self.dist.set_disttype('d_star')
        matrix = distmatrix.create(self.pep_records.id_list, self.dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.2250000 0.3750000 0.5000000",
            "seq2       0.2250000 0.0000000 0.2631579 0.4736842",
            "seq3       0.3750000 0.2631579 0.0000000 0.4857143",
            "seq4       0.5000000 0.4736842 0.4857143 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance_d1(self):
        self.dist.set_disttype('d1')
        matrix = distmatrix.create(self.pep_records.id_list, self.dist)
        exp = [
            "   4",
            "seq1       0 16 28 23",
            "seq2       16 0 19 21",
            "seq3       28 19 0 19",
            "seq4       23 21 19 0"
        ]
        self.assertEqual(matrix.format(0), "\n".join(exp))

    def test_distance_d1_star(self):
        self.dist.set_disttype('d1_star')
        matrix = distmatrix.create(self.pep_records.id_list, self.dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.3404255 0.5283019 0.5348837",
            "seq2       0.3404255 0.0000000 0.4042553 0.5121951",
            "seq3       0.5283019 0.4042553 0.0000000 0.5135135",
            "seq4       0.5348837 0.5121951 0.5135135 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance_d1_star2(self):
        self.dist.set_disttype('d1_star2')
        matrix = distmatrix.create(self.pep_records.id_list, self.dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.3404255 0.5436893 0.5609756",
            "seq2       0.3404255 0.0000000 0.4130435 0.5384615",
            "seq3       0.5436893 0.4130435 0.0000000 0.5205479",
            "seq4       0.5609756 0.5384615 0.5205479 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_set_disttype_throws_exception(self):
        with self.assertRaises(Exception) as context:
            self.dist.set_disttype('nonexitent')
        self.assertIn('unknown disttype', str(context.exception))


if __name__ == '__main__':
    unittest.main()
