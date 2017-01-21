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

    def test_distance_yule(self):
        dist = word_bool_distance.Distance(self.vector, 'yule')
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.3213773 1.7495396 0.1451613",
            "seq2       0.3213773 0.0000000 0.4636488 0.0955414",
            "seq3       1.7495396 0.4636488 0.0000000 0.0000000",
            "seq4       0.1451613 0.0955414 0.0000000 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance_rogerstanimoto(self):
        dist = word_bool_distance.Distance(self.vector, 'rogerstanimoto')
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.4313725 0.7096774 0.6324786",
            "seq2       0.4313725 0.0000000 0.4905660 0.5585586",
            "seq3       0.7096774 0.4905660 0.0000000 0.5321101",
            "seq4       0.6324786 0.5585586 0.5321101 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance_russellrao(self):
        dist = word_bool_distance.Distance(self.vector, 'russellrao')
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.4375000 0.5750000 0.7125000",
            "seq2       0.4375000 0.0000000 0.5000000 0.7125000",
            "seq3       0.5750000 0.5000000 0.0000000 0.7000000",
            "seq4       0.7125000 0.7125000 0.7000000 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance_sokalsneath(self):
        dist = word_bool_distance.Distance(self.vector, 'sokalsneath')
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.4943820 0.7213115 0.7628866",
            "seq2       0.4943820 0.0000000 0.5652174 0.7294118",
            "seq3       0.7213115 0.5652174 0.0000000 0.7073171",
            "seq4       0.7628866 0.7294118 0.7073171 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance_sokalsneath_throws_exception(self):
        vector = np.array([[False, False, False], [False, False, False]])
        dist = word_bool_distance.Distance(vector, 'sokalsneath')
        with self.assertRaises(ValueError) as context:
            matrix = distmatrix.create(['s1', 's2'], dist)
        exp = 'vectors that are entirely false'
        self.assertIn(exp, str(context.exception))

    def test_distance_sokalmichener(self):
        dist = word_bool_distance.Distance(self.vector, 'sokalmichener')
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.4313725 0.7096774 0.6324786",
            "seq2       0.4313725 0.0000000 0.4905660 0.5585586",
            "seq3       0.7096774 0.4905660 0.0000000 0.5321101",
            "seq4       0.6324786 0.5585586 0.5321101 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance_kulsinski(self):
        dist = word_bool_distance.Distance(self.vector, 'kulsinski')
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.5588235 0.7258065 0.8034188",
            "seq2       0.5588235 0.0000000 0.6226415 0.7927928",
            "seq3       0.7258065 0.6226415 0.0000000 0.7798165",
            "seq4       0.8034188 0.7927928 0.7798165 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))


if __name__ == '__main__':
    unittest.main()
