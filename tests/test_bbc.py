import unittest

from alfpy import bbc
from alfpy.utils import distmatrix

from . import utils


class VectorTest(unittest.TestCase, utils.ModulesCommonTest):
    """Shared methods and tests for creating BBC vectors."""

    def __init__(self, *args, **kwargs):
        super(VectorTest, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()

    def test_base_base_correlation_on_dna_k1(self):
        seq = 'CTAGGGAACATACCA'
        vec = bbc.base_base_correlation(seq, 1, utils.ALPHABET_DNA)
        exp = [[-0.09531944216197286, 0.08354249482243745,
                -0.011648745118694886, -0.03537936275852141,
                0.08354249482243745, -0.07114045890800086,
                -0.10696198320987654, 0.07735938746326047,
                0.06569448291953908, -0.10696198320987654,
                0.14092829842677065, -0.05337075358024692,
                0.036771706024547365, 0.0058396884306785996,
                0.02369886445107525, -0.0355666929797287]]
        self.assertEqual(vec.shape, (1, 16))
        self.assertEqual(vec.tolist(), exp)

    def test_base_base_correlation_on_dna_k1_alphabet_None(self):
        seq = 'CTAGGGAACATACCA'
        vec = bbc.base_base_correlation(seq, 1)
        exp = [[-0.09531944216197286, 0.08354249482243745,
                -0.011648745118694886, -0.03537936275852141,
                0.08354249482243745, -0.07114045890800086,
                -0.10696198320987654, 0.07735938746326047,
                0.06569448291953908, -0.10696198320987654,
                0.14092829842677065, -0.05337075358024692,
                0.036771706024547365, 0.0058396884306785996,
                0.02369886445107525, -0.0355666929797287]]
        self.assertEqual(vec.shape, (1, 16))
        self.assertEqual(vec.tolist(), exp)

    def test_base_base_correlation_on_ambiguous_dna_k1(self):
        """Base_base_correlation function works when a sequence contains
        characters (e.g. Ns) not included in alphabet.

        """
        dna = 'CTAGGGAACATACCANNNNNNNN'
        vec = bbc.base_base_correlation(dna, 1, utils.ALPHABET_DNA)
        exp = [[-0.09531944216197286, 0.08354249482243745,
                -0.011648745118694886, -0.03537936275852141,
                0.08354249482243745, -0.07114045890800086,
                -0.10696198320987654, 0.07735938746326047,
                0.06569448291953908, -0.10696198320987654,
                0.14092829842677065, -0.05337075358024692,
                0.036771706024547365, 0.0058396884306785996,
                0.02369886445107525, -0.0355666929797287]]
        self.assertEqual(vec.shape, (1, 16))
        self.assertEqual(vec.tolist(), exp)

    def test_base_base_correlation_on_dna_k2(self):
        vec = bbc.base_base_correlation('ATGCATGC', 2, utils.ALPHABET_DNA)
        exp = [[-0.18820953369140625, 0.2506928100585938,
                0.16536938702618637, 0.1089601433311885,
                -0.044584379946872726, -0.18820953369140625,
                0.014663513183593768, -0.020154882360387746,
                -0.020154882360387746, 0.1089601433311885,
                -0.18820953369140625, 0.014663513183593768,
                0.014663513183593768, 0.16536938702618637,
                0.1089601433311885, -0.18820953369140625]]
        self.assertEqual(vec.shape, (1, 16))
        self.assertEqual(vec.tolist(), exp)

    def test_base_base_correlation_on_pep_k1(self):
        """Base_base_correlation function also handles protein sequences."""

        seq = 'MFTDNAKIIIVQLNASVEINCTRPNNNTR'
        vec = bbc.base_base_correlation(seq, 1, utils.ALPHABET_PEP)
        md5 = utils.calc_md5(vec)
        self.assertEqual(vec.shape, (1, 400))
        self.assertEqual(md5, 'fbcf029bcc745b076427ebfb9562ae49')

    def test_base_base_correlation_throws_exception(self):
        # Function raises an exception if the sequence
        # is too short (i.e. len(s) - 2 < k)
        with self.assertRaises(Exception) as context:
            bbc.base_base_correlation('ACT', 2, utils.ALPHABET_DNA)
        self.assertIn('Sequence too short', str(context.exception))

    def test_create_vectors_on_dna_k1(self):
        vec = bbc.create_vectors(self.dna_records, 1, utils.ALPHABET_DNA)
        md5 = utils.calc_md5(vec)
        self.assertEqual(vec.shape, (3, 16))
        self.assertEqual(md5, 'e99bc40356b00e04fd858a665af597ec')

    def test_create_vectors_on_pep_k1(self):
        vec = bbc.create_vectors(self.pep_records, 1, utils.ALPHABET_PEP)
        md5 = utils.calc_md5(vec)
        self.assertEqual(vec.shape, (4, 400))
        self.assertEqual(md5, 'd901d93d0d71102d0727633a67fc14a3')

    def test_create_vectors_throws_exception(self):
        with self.assertRaises(Exception) as context:
            bbc.create_vectors(self.dna_records, 60, utils.ALPHABET_DNA)
        self.assertIn('Sequence too short', str(context.exception))


class DistanceTest(unittest.TestCase, utils.ModulesCommonTest):
    """Shared methods and tests for Distances calculations."""

    def __init__(self, *args, **kwargs):
        super(DistanceTest, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()
        self.vector = bbc.create_vectors(self.dna_records, 10,
                                         utils.ALPHABET_DNA)

    def test_distance_dna_euclidnorm(self):
        dist = bbc.Distance(self.vector)
        matrix = distmatrix.create(self.dna_records.id_list, dist)
        exp = [
            "   3",
            "seq1       0.0000000 1.0227476 1.9351116",
            "seq2       1.0227476 0.0000000 1.4469591",
            "seq3       1.9351116 1.4469591 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance_dna_google(self):
        dist = bbc.Distance(self.vector, 'google')
        matrix = distmatrix.create(self.dna_records.id_list, dist)
        exp = [
            "   3",
            "seq1       0.0000000 73.1311144 37.1219467",
            "seq2       73.1311144 0.0000000 33.2221873",
            "seq3       37.1219467 33.2221873 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))

    def test_distance_pep_google(self):
        vector = bbc.create_vectors(self.pep_records, 10, utils.ALPHABET_PEP)
        dist = bbc.Distance(vector, 'google')
        matrix = distmatrix.create(self.pep_records.id_list, dist)
        exp = [
            "   4",
            "seq1       0.0000000 862.9400129 912.1886801 233.2058119",
            "seq2       862.9400129 0.0000000 694.9729647 237.9720958",
            "seq3       912.1886801 694.9729647 0.0000000 225.2262758",
            "seq4       233.2058119 237.9720958 225.2262758 0.0000000"
        ]
        self.assertEqual(matrix.format(), "\n".join(exp))


if __name__ == '__main__':
    unittest.main()
