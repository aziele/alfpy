import numpy as np
import hashlib
import os
import unittest

from alfpy import graph2Ddna
from alfpy.utils import distmatrix

from . import utils


class VectorTest(unittest.TestCase, utils.ModulesCommonTest):
    """Shared methods and tests for creating Graph 2D DNA vectors."""

    def __init__(self, *args, **kwargs):
        super(VectorTest, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()

    def test_2DSGraphVector(self):
        seq = 'CTAGGGAACATACCA'
        vec = graph2Ddna._2DSGraphVector(seq)

        exp = [2.99197183, -8.04298066, 9.16666667, -5.78272208,
               6.5, -1.75064326, 5, -2.92241364, 9.25, -3.81343559]
        self.assertTrue(np.allclose(vec, np.array(exp)))

    def test_2DMGraphVector_ndim10(self):
        seq = 'CTAGGGAACATACCA'
        vec = graph2Ddna._2DMGraphVector(seq, 10)
        exp = [15, 12.14790682, 13.5804606, 15.88980624, 19.16010756,
               23.55763468, 29.38627489, 37.08035601, 47.23633868,
               60.66394053]
        self.assertEqual(vec.shape, (10,))
        self.assertTrue(np.allclose(vec, np.array(exp)))

    def test_2DMGraphVector_ndim5(self):
        seq = 'CTAGGGAACATACCA'
        vec = graph2Ddna._2DMGraphVector(seq, 5)
        exp = [15, 12.14790682, 13.5804606, 15.88980624, 19.16010756]
        self.assertEqual(vec.shape, (5,))
        self.assertTrue(np.allclose(vec, np.array(exp)))

    def test_2DNGraphVector(self):
        seq = 'CTAGGGAACATACCA'
        vec = graph2Ddna._2DNGraphVector(seq)
        md5 = utils.calc_md5(vec)
        self.assertEqual(len(vec), 48)
        self.assertEqual(md5, '44829cc0277531646d656cdaacd3ae94')

    def test_create_2DSGraphVectors(self):
        data = graph2Ddna.create_2DSGraphVectors(self.dna_records)
        md5 = utils.calc_md5(data)
        self.assertEqual(md5, 'e2399897bb7eaa5ca3a81c84e2eeac84')

    def test_create_2DMGraphVectors(self):
        data = graph2Ddna.create_2DMGraphVectors(self.dna_records, 10)
        md5 = utils.calc_md5(data)
        self.assertEqual(md5, '8c7d4dca912aeaf7c88d325799dadf00')

    def test_create_2DNGraphVectors(self):
        data = graph2Ddna.create_2DNGraphVectors(self.dna_records)
        md5 = utils.calc_md5(data)
        self.assertEqual(md5, 'd24c8508fb6f2b833cfe45f7af38be34')


class DistanceTest(unittest.TestCase, utils.ModulesCommonTest):
    """Shared methods and tests for Distances calculations."""

    def __init__(self, *args, **kwargs):
        super(DistanceTest, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()

    def test_distance_2DSG(self):
        data = graph2Ddna.create_2DSGraphVectors(self.dna_records)
        dist = graph2Ddna.Distance(data)
        matrix = distmatrix.create(self.dna_records.id_list, dist)
        exp = [
            '   3',
            'seq1       0.0000000 9.4762599 14.6585286',
            'seq2       9.4762599 0.0000000 6.7199568',
            'seq3       14.6585286 6.7199568 0.0000000',
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance_2DMG(self):
        data = graph2Ddna.create_2DMGraphVectors(self.dna_records, 10)
        dist = graph2Ddna.Distance(data)
        matrix = distmatrix.create(self.dna_records.id_list, dist)
        exp = [
            '   3',
            'seq1       0.0000000 22.2449494 55.9753388',
            'seq2       22.2449494 0.0000000 34.2064423',
            'seq3       55.9753388 34.2064423 0.0000000'
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance_2DNG(self):
        data = graph2Ddna.create_2DNGraphVectors(self.dna_records)
        dist = graph2Ddna.Distance(data)
        matrix = distmatrix.create(self.dna_records.id_list, dist)
        exp = [
            '   3',
            'seq1       0.0000000 10.3711467 15.1355787',
            'seq2       10.3711467 0.0000000 7.8973545',
            'seq3       15.1355787 7.8973545 0.0000000'
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))


if __name__ == '__main__':
    unittest.main()
