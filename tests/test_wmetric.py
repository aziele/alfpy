import unittest

from alfpy import wmetric
from alfpy.utils import distmatrix
from alfpy.utils.data import subsmat

from . import utils


class VectorTest(unittest.TestCase):

    def test_count_seq_chars(self):
        seq = 'MKSTGWHFSG'
        l = wmetric.count_seq_chars(seq, utils.ALPHABET_PEP)
        expl = [0, 0, 0, 0, 1, 2, 1, 0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0]
        self.assertEqual(l, expl)

    def test_count_seq_chars_pep_ambiguous(self):
        seq = 'MKSTGWXXXXXXXOOOOOOOHFSG'
        l = wmetric.count_seq_chars(seq, utils.ALPHABET_PEP)
        expl = [0, 0, 0, 0, 1, 2, 1, 0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0]
        self.assertEqual(l, expl)

    def test_freq_seq_chars(self):
        seq = 'MKSTGWXXXXXXXOOOOOOOHFSG'
        l = wmetric.count_seq_chars(seq, utils.ALPHABET_PEP)
        freq = wmetric.freq_seq_chars(l)
        expfreq = [0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.1, 0.0, 0.1, 0.0,
                   0.1, 0.0, 0.0, 0.0, 0.2, 0.1, 0.0, 0.1, 0.0, 0.0]
        self.assertEqual(freq, expfreq)


class DistanceTest(unittest.TestCase, utils.ModulesCommonTest):

    def __init__(self, *args, **kwargs):
        super(DistanceTest, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()

    def test_wmetric_blosum62(self):
        # The result of this method is identical to that from decaf+py.
        matrix = subsmat.get('blosum62')
        dist = wmetric.Distance(self.pep_records, matrix)
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        data = ['   4',
                'seq1       0.0000000 0.0392559 0.0783026 0.1261381',
                'seq2       0.0392559 0.0000000 0.0377364 0.1166475',
                'seq3       0.0783026 0.0377364 0.0000000 0.1677386',
                'seq4       0.1261381 0.1166475 0.1677386 0.0000000']
        self.assertEqual(matrix.format(), "\n".join(data))

    def test_wmetric_pam250(self):
        matrix = subsmat.get('pam250')
        dist = wmetric.Distance(self.pep_records, matrix)
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        data = ['   4',
                'seq1       0.0000000 0.0289700 0.0467580 0.0353781',
                'seq2       0.0289700 0.0000000 0.0227122 0.0372699',
                'seq3       0.0467580 0.0227122 0.0000000 0.0578383',
                'seq4       0.0353781 0.0372699 0.0578383 0.0000000']
        self.assertEqual(matrix.format(), "\n".join(data))


if __name__ == '__main__':
    unittest.main()
