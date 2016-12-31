import unittest

from alfpy import word_d2
from alfpy import word_pattern
from alfpy import word_vector
from alfpy.utils import distmatrix

from . import utils


class DistanceTest(unittest.TestCase, utils.ModulesCommonTest):

    def __init__(self, *args, **kwargs):
        super(DistanceTest, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()
        self.patterns = []
        self.counts = []
        self.freqs = []
        for i in range(1, 5):
            p = word_pattern.create(self.pep_records.seq_list, i)
            self.patterns.append(p)
            c = word_vector.Counts(self.pep_records.length_list, p)
            self.counts.append(c)
            f = word_vector.Freqs(self.pep_records.length_list, p)
            self.freqs.append(f)

    def test_counts_from1_to4(self):
        dist = word_d2.Distance(self.counts)
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            '   4',
            'seq1       0 130 236 286',
            'seq2       130 0 142 258',
            'seq3       236 142 0 212',
            'seq4       286 258 212 0'
        ]
        self.assertEqual(matrix.format(decimal_places=0), "\n".join(exp))

    def test_freqs_from1_to4(self):
        dist = word_d2.Distance(self.freqs)
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            '   4',
            'seq1       0.0000000 0.0313590 0.0573154 0.1020235',
            'seq2       0.0313590 0.0000000 0.0373677 0.0907196',
            'seq3       0.0573154 0.0373677 0.0000000 0.0870581',
            'seq4       0.1020235 0.0907196 0.0870581 0.0000000'

        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_counts_from1_to1(self):
        dist = word_d2.Distance([self.counts[0]])
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            '   4',
            'seq1       0 37 57 140',
            'seq2       37 0 28 137',
            'seq3       57 28 0 111',
            'seq4       140 137 111 0'
        ]
        self.assertEqual(matrix.format(decimal_places=0), "\n".join(exp))

    def test_freqs_from1_to4_d2_squareroot(self):
        dist = word_d2.Distance(self.freqs)
        dist.set_disttype('d2_squareroot')
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0.0000000 0.1770847 0.2394063 0.3194113",
            "seq2       0.1770847 0.0000000 0.1933073 0.3011969",
            "seq3       0.2394063 0.1933073 0.0000000 0.2950560",
            "seq4       0.3194113 0.3011969 0.2950560 0.0000000"

        ]
        self.assertEqual(matrix.format(), "\n".join(exp))


if __name__ == '__main__':
    unittest.main()
