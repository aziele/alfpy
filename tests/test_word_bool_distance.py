import numpy as np
import unittest

from alfpy import word_pattern
from alfpy import word_vector
from alfpy import word_bool_distance
from alfpy.utils import distmatrix

from . import utils


class Test(unittest.TestCase, utils.ModulesCommonTest):

    def __init__(self, *args, **kwargs):
        super(Test, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()
        self.p = word_pattern.create(self.pep_records.seq_list, 2)
        self.vector = word_vector.Bools(self.pep_records.length_list, self.p)

    def test_nbool_correspond_ft_tf(self):
        u = np.array([True, False, True])
        v = np.array([True, True, False])
        vec = word_bool_distance._nbool_correspond_ft_tf(u, v)
        self.assertEqual(vec, (1, 1))

    def test_nbool_correspond_ft_tf_u_equals_v(self):
        u = np.array([True, True, False])
        v = np.array([True, True, False])
        vec = word_bool_distance._nbool_correspond_ft_tf(u, v)
        self.assertEqual(vec, (0, 0))

    def test_nbool_correspond_all(self):
        u = np.array([True, False, True])
        v = np.array([True, True, False])
        vec = word_bool_distance._nbool_correspond_all(u, v)
        self.assertEqual(vec, (0, 1, 1, 1))

    def test_distance_dice(self):
        dist = word_bool_distance.Distance(self.vector, 'dice')
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.1964286 0.3928571 0.4457831",
            "seq2       0.1964286 0.0000000 0.2452830 0.4025974",
            "seq3       0.3928571 0.2452830 0.0000000 0.3766234",
            "seq4       0.4457831 0.4025974 0.3766234 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance_hamming(self):
        dist = word_bool_distance.Distance(self.vector, 'hamming')
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.2750000 0.5500000 0.4625000",
            "seq2       0.2750000 0.0000000 0.3250000 0.3875000",
            "seq3       0.5500000 0.3250000 0.0000000 0.3625000",
            "seq4       0.4625000 0.3875000 0.3625000 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance_jaccard(self):
        dist = word_bool_distance.Distance(self.vector, 'jaccard')
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.3283582 0.5641026 0.6166667",
            "seq2       0.3283582 0.0000000 0.3939394 0.5740741",
            "seq3       0.5641026 0.3939394 0.0000000 0.5471698",
            "seq4       0.6166667 0.5740741 0.5471698 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))


if __name__ == '__main__':
    unittest.main()
