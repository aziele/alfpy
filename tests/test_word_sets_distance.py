import numpy as np
import unittest

from alfpy import word_pattern
from alfpy import word_sets_distance
from alfpy.utils import distmatrix

from . import utils


class Test(unittest.TestCase, utils.ModulesCommonTest):

    def __init__(self, *args, **kwargs):
        super(Test, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()
        self.p = word_pattern.create(self.pep_records.seq_list, 2)

    def test_getwords(self):
        words = word_sets_distance._getwords('ATGCGTA', 2)
        self.assertSetEqual(words, set(['GT', 'CG', 'GC', 'AT', 'TG', 'TA']))

    def test_distance_dice(self):
        # The result of this function is identical
        # to the Dice distance implemented in word_bool_distance.
        dist = word_sets_distance.Distance(self.pep_records, 2, 'dice')
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
        dist = word_sets_distance.Distance(self.pep_records, 2, 'hamming')
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0 22 44 37",
            "seq2       22 0 26 31",
            "seq3       44 26 0 29",
            "seq4       37 31 29 0"
        ]
        self.assertEqual(matrix.format(0), "\n".join(exp))

    def test_distance_jaccard(self):
        # The result of this function is identical
        # to the Jaccard distance implemented in word_bool_distance.
        dist = word_sets_distance.Distance(self.pep_records, 2, 'jaccard')
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
