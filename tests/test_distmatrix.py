import numpy as np
import os
import unittest

from alfpy import word_distance
from alfpy.utils import distmatrix

from . import utils


class TestDistMatrix(unittest.TestCase):

    def setUp(self):
        id_list = ['seq1', 'seq2', 'seq3']
        data = np.array([[0, 0.3531587, 0.35509333],
                         [0.3531587, 0, 0.295394],
                         [0.35509333, 0.295394, 0.]
                         ])
        self.matrix = distmatrix.Matrix(id_list, data)
        self.output_filename = utils.get_test_data('distmatrix.txt')

    def test_format(self):
        exp = [
            '   3',
            'seq1       0.0000000 0.3531587 0.3550933',
            'seq2       0.3531587 0.0000000 0.2953940',
            'seq3       0.3550933 0.2953940 0.0000000'
        ]
        self.assertEqual(self.matrix.format(), "\n".join(exp))

    def test_format_decimal3(self):
        exp = [
            '   3',
            'seq1       0.000 0.353 0.355',
            'seq2       0.353 0.000 0.295',
            'seq3       0.355 0.295 0.000'
        ]
        self.assertEqual(self.matrix.format(3), "\n".join(exp))

    def test_min(self):
        self.assertEqual(self.matrix.min(), 0)

    def test_max(self):
        self.assertEqual(self.matrix.max(), 0.35509332999999998)

    def test_is_zero(self):
        self.assertFalse(self.matrix.is_zero())

    def test_normalize(self):
        self.matrix.normalize()
        exp = [
            "   3",
            "seq1       0.0000000 0.9945518 1.0000000",
            "seq2       0.9945518 0.0000000 0.8318771",
            "seq3       1.0000000 0.8318771 0.0000000",
        ]
        self.assertEqual(self.matrix.format(), "\n".join(exp))

    def test_write_to_file_phylip(self):
        oh = open(self.output_filename, 'w')
        self.matrix.write_to_file(oh)
        oh.close()
        fh = open(self.output_filename)
        result = fh.read()
        fh.close()
        os.remove(self.output_filename)
        exp = [
            '   3',
            'seq1       0.0000000 0.3531587 0.3550933',
            'seq2       0.3531587 0.0000000 0.2953940',
            'seq3       0.3550933 0.2953940 0.0000000\n'
        ]
        self.assertEqual(result, "\n".join(exp))

    def test_write_to_file_pairwise(self):
        oh = open(self.output_filename, 'w')
        self.matrix.write_to_file(oh, 'pairwise')
        oh.close()
        fh = open(self.output_filename)
        result = fh.read()
        fh.close()
        os.remove(self.output_filename)
        exp = [
            "seq1\tseq2\t0.3531587",
            "seq1\tseq3\t0.3550933",
            "seq2\tseq3\t0.2953940\n"
        ]
        self.assertEqual(result, "\n".join(exp))

    def test_write_to_file_pairwise_decimal3(self):
        oh = open(self.output_filename, 'w')
        self.matrix.write_to_file(oh, 'pairwise', 3)
        oh.close()
        fh = open(self.output_filename)
        result = fh.read()
        fh.close()
        os.remove(self.output_filename)
        exp = [
            "seq1\tseq2\t0.353",
            "seq1\tseq3\t0.355",
            "seq2\tseq3\t0.295\n"
        ]
        self.assertEqual(result, "\n".join(exp))

    def test_iter(self):
        exp = [(0, 1, 'seq1', 'seq2', 0.35315869999999999),
               (0, 2, 'seq1', 'seq3', 0.35509332999999998),
               (1, 2, 'seq2', 'seq3', 0.29539399999999999)]
        self.assertEqual(list(self.matrix), exp)

    def test_create_matrix(self):
        l = [[3, 6, 4, 1, 3, 4, 3, 0, 1, 1, 6, 4, 5, 0, 3, 4],
             [0, 3, 0, 3, 0, 0, 0, 2, 9, 0, 3, 3, 0, 6, 3, 6],
             [9, 0, 0, 3, 0, 0, 0, 2, 6, 0, 3, 3, 0, 3, 3, 3]]
        vector = np.array(l)
        dist = word_distance.Distance(vector, 'minkowski')
        id_list = ['seq1', 'seq2', 'seq3']
        matrix = distmatrix.create(id_list, dist)
        exp = [
            '   3',
            'seq1       0.0000000 14.6969385 14.1774469',
            'seq2       14.6969385 0.0000000 10.8166538',
            'seq3       14.1774469 10.8166538 0.0000000'
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_highcharts(self):
        self.assertEqual(len(self.matrix.highcharts()), 3)

    def test_read_highcharts_matrix(self):
        id_list = ['seq1', 'seq2', 'seq3']
        data = [[0, 1, 0.35, 0.19], [0, 2, 1.0, 0.55], [1, 2, 0.88, 0.48]]
        matrix = distmatrix.read_highcharts_matrix(id_list, data)
        md5 = utils.calc_md5(matrix.format())
        self.assertEqual(md5, "476c8f5d284a84ee3c7c419bde2d7658")
        

if __name__ == '__main__':
    unittest.main()
