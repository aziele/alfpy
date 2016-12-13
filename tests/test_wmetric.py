import unittest

from alfpy import wmetric
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords
from alfpy.utils.data import subsmat


class TestWmetric(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestWmetric, self).__init__(*args, **kwargs)
        fh = open('data/pep.fa')
        self.seq_records = seqrecords.read_fasta(fh)
        fh.close()
        self.alphabet = 'ACDEFGHIKLMNPQRSTVWY'

    def test_count_seq_chars(self):
        seq = 'MKSTGWXXXXXXXOOOOOOOHFSG'
        l = wmetric.count_seq_chars(seq, self.alphabet)
        freq = wmetric.freq_seq_chars(l)
        expfreq = [0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.1, 0.0, 0.1, 0.0,
                   0.1, 0.0, 0.0, 0.0, 0.0, 0.2, 0.1, 0.0, 0.1, 0.0]
        self.assertEqual(freq, expfreq)

    def test_freq_seq_chars(self):
        seq = 'MKSTGWXXXXXXXOOOOOOOHFSG'
        l = wmetric.count_seq_chars(seq, self.alphabet)
        expl = [0, 0, 0, 0, 1, 2, 1, 0, 1, 0, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0]
        self.assertEqual(l, expl)

    def test_wmetric_blosum62(self):
        # agrees with decaf+py
        matrix = subsmat.get('blosum62')
        dist = wmetric.Distance(self.seq_records, matrix)
        matrix = distmatrix.create(self.seq_records.id_list, dist)
        data = ['   4',
                'seq1       0.0000000 0.0392559 0.0783026 0.1261381',
                'seq2       0.0392559 0.0000000 0.0377364 0.1166475',
                'seq3       0.0783026 0.0377364 0.0000000 0.1677386',
                'seq4       0.1261381 0.1166475 0.1677386 0.0000000\n']
        self.assertEqual(matrix.format(), "\n".join(data))

    def test_wmetric_pam250(self):
        matrix = subsmat.get('pam250')
        dist = wmetric.Distance(self.seq_records, matrix)
        matrix = distmatrix.create(self.seq_records.id_list, dist)
        data = ['   4',
                'seq1       0.0000000 0.0289700 0.0467580 0.0353781',
                'seq2       0.0289700 0.0000000 0.0227122 0.0372699',
                'seq3       0.0467580 0.0227122 0.0000000 0.0578383',
                'seq4       0.0353781 0.0372699 0.0578383 0.0000000\n']
        self.assertEqual(matrix.format(), "\n".join(data))


if __name__ == '__main__':
    unittest.main()
