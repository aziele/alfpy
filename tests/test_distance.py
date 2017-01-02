import unittest

from alfpy import word_pattern
from alfpy import word_vector
from alfpy.utils import distance
from alfpy.utils import distmatrix

from . import utils


class DistanceTest(unittest.TestCase, utils.ModulesCommonTest):

    def __init__(self, *args, **kwargs):
        super(DistanceTest, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()
        self.pattern = word_pattern.create(self.dna_records.seq_list, 2)
        self.counts = word_vector.Counts(self.dna_records.length_list,
                                         self.pattern)
        self.freqs = word_vector.Freqs(self.dna_records.length_list,
                                       self.pattern)

    def test_euclid_squared_counts(self):
        # The result of this method is identical to that from decaf+py.
        dist = distance.Distance(self.counts, 'euclid_squared')
        matrix = distmatrix.create(self.dna_records.id_list, dist)
        data = ['   3',
                'seq1       0.0000000 57.0000000 30.0000000',
                'seq2       57.0000000 0.0000000 19.0000000',
                'seq3       30.0000000 19.0000000 0.0000000']
        self.assertEqual(matrix.format(), "\n".join(data))

    def test_euclid_squared_freqs(self):
        # The result of this method is identical to that from decaf+py.
        dist = distance.Distance(self.freqs, 'euclid_squared')
        matrix = distmatrix.create(self.dna_records.id_list, dist)
        data = ['   3',
                'seq1       0.0000000 0.1416402 0.0641298',
                'seq2       0.1416402 0.0000000 0.0677565',
                'seq3       0.0641298 0.0677565 0.0000000']
        self.assertEqual(matrix.format(), "\n".join(data))

    def test_euclid_norm_counts(self):
        # The result of this method is identical to that from decaf+py.
        dist = distance.Distance(self.counts, 'euclid_norm')
        matrix = distmatrix.create(self.dna_records.id_list, dist)
        data = ['   3',
                'seq1       0.0000000 7.5498344 5.4772256',
                'seq2       7.5498344 0.0000000 4.3588989',
                'seq3       5.4772256 4.3588989 0.0000000']
        self.assertEqual(matrix.format(), "\n".join(data))

    def test_euclid_norm_freqs(self):
        # The result of this method is identical to that from decaf+py.
        dist = distance.Distance(self.freqs, 'euclid_norm')
        matrix = distmatrix.create(self.dna_records.id_list, dist)
        data = ['   3',
                'seq1       0.0000000 0.3763512 0.2532387',
                'seq2       0.3763512 0.0000000 0.2603008',
                'seq3       0.2532387 0.2603008 0.0000000']
        self.assertEqual(matrix.format(), "\n".join(data))

    def test_google_freqs(self):
        dist = distance.Distance(self.freqs, 'google')
        matrix = distmatrix.create(self.dna_records.id_list, dist)
        data = ['   3',
                'seq1       0.0000000 0.6078431 0.3809524',
                'seq2       0.6078431 0.0000000 0.3949580',
                'seq3       0.3809524 0.3949580 0.0000000']
        self.assertEqual(matrix.format(), "\n".join(data))

    def test_get_disttypes(self):
        distlist = distance.Distance.get_disttypes()
        exp = ['euclid_norm', 'euclid_squared', 'google']
        self.assertListEqual(distlist, exp)

    def test_set_disttypes_throws_exception(self):
        dist = distance.Distance(self.freqs, 'google')
        with self.assertRaises(Exception) as context:
            dist.set_disttype('nonexistent')
        self.assertIn('unknown disttype', str(context.exception))

if __name__ == '__main__':
    unittest.main()
