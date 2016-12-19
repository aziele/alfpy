import os
import unittest

from alfpy import fcgr
from alfpy.utils import distmatrix

from . import utils


class VectorTest(unittest.TestCase, utils.ModulesCommonTest):
    """Shared methods and tests for creating BBC vectors."""

    def __init__(self, *args, **kwargs):
        super(VectorTest, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()

    def test_fcgr_vector1(self):
        vec = fcgr.fcgr_vector('CTAGGGAACATACCA', 1)
        self.assertEqual(vec, [3.0, 6.0, 3.0])

    def test_fcgr_vector2(self):
        vec = fcgr.fcgr_vector('CTAGGGAACATACCA', 2)
        exp = [0.0, 0.0, 2.0, 2.0, 1.0, 1.0, 0.0, 2.0,
               1.0, 2.0, 0.0, 1.0, 2.0, 1.0, 0.0]
        self.assertEqual(vec, exp)

    def test_fcgr_vector3(self):
        vec = fcgr.fcgr_vector('CTAGGGAACATACCXXA', 1)
        self.assertEqual(vec, [3.0, 6.0, 3.0])

    def test_create_vectors(self):
        vecs = fcgr.create_vectors(self.dna_records, 2)
        exp = [[0, 3, 1, 4, 0, 1, 1, 1, 1, 1, 3, 2, 4, 1, 1],
               [0, 0, 4, 1, 2, 2, 0, 0, 1, 4, 0, 0, 3, 1, 1],
               [0, 0, 2, 2, 1, 1, 0, 2, 1, 2, 0, 1, 2, 1, 0]]
        self.assertEqual(vecs.tolist(), exp)


class DistanceTest(unittest.TestCase, utils.ModulesCommonTest):
    """Shared methods and tests for Distances calculations."""

    def __init__(self, *args, **kwargs):
        super(DistanceTest, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()

    def test_distance1(self):
        vecs = fcgr.create_vectors(self.dna_records, 2)
        dist = fcgr.Distance(vecs)
        matrix = distmatrix.create(self.dna_records.id_list, dist)
        exp = [
            "   3",
            "seq1       0.0000000 7.5498344 5.7445626",
            "seq2       7.5498344 0.0000000 4.2426407",
            "seq3       5.7445626 4.2426407 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance2(self):
        vecs = fcgr.create_vectors(self.dna_records, 2)
        dist = fcgr.Distance(vecs, 'google')
        matrix = distmatrix.create(self.dna_records.id_list, dist)
        exp = [
            "   3",
            "seq1       0.0000000 0.5833333 0.5416667",
            "seq2       0.5833333 0.0000000 0.4210526",
            "seq3       0.5416667 0.4210526 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))


if __name__ == '__main__':
    unittest.main()
