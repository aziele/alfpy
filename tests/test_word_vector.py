import unittest

from alfpy import word_pattern
from alfpy import word_vector

from . import utils


class WordVectorTest(unittest.TestCase, utils.ModulesCommonTest):

    def __init__(self, *args, **kwargs):
        super(WordVectorTest, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()
        self.pattern1 = word_pattern.create(self.dna_records.seq_list, 1)
        self.pattern2 = word_pattern.create(self.dna_records.seq_list, 2)

    def test_counts_pattern1(self):
        counts = word_vector.Counts(self.dna_records.length_list,
                                    self.pattern1)
        exp = ["A\t8 4 6",
               "G\t6 6 3",
               "C\t6 3 4",
               "T\t5 5 2"]
        lengths = [25, 18, 15]
        self.assertEqual(counts.format(decimal_places=0), "\n".join(exp))
        # Counts in a sequence should sum to sequence length.
        for i in range(len(counts.data)):
            self.assertEqual(sum(counts.data[i]), lengths[i])

    def test_counts_pattern2(self):
        counts = word_vector.Counts(self.dna_records.length_list,
                                    self.pattern2)
        exp = [
            "TA\t3 3 2",
            "AC\t4 1 2",
            "GG\t1 4 2",
            "AG\t1 2 1",
            "CT\t0 3 1",
            "AA\t2 0 1",
            "AT\t1 1 1",
            "CA\t1 0 2",
            "CC\t2 0 1",
            "CG\t3 0 0",
            "GA\t1 1 1",
            "GT\t3 0 0",
            "TT\t1 1 0",
            "TC\t0 1 0",
            "TG\t1 0 0"
        ]
        self.assertEqual(counts.format(decimal_places=0), "\n".join(exp))
        for i in range(len(counts.data)):
            self.assertEqual(sum(counts.data[i]), counts.seq_lengths[i] - 1)

    def test_freqs_pattern1(self):
        freqs = word_vector.Freqs(self.dna_records.length_list,
                                  self.pattern1)

        exp = [
            "A\t0.320 0.222 0.400",
            "G\t0.240 0.333 0.200",
            "C\t0.240 0.167 0.267",
            "T\t0.200 0.278 0.133",
        ]
        self.assertEqual(freqs.format(), "\n".join(exp))
        # Freqs in a given sequence should sum to 1.
        for seqrow in freqs.data:
            self.assertEqual('{:.3f}'.format(sum(seqrow)), '1.000')

    def test_freqs_pattern2(self):
        freqs = word_vector.Freqs(self.dna_records.length_list,
                                  self.pattern2)
        for seqrow in freqs.data:
            self.assertEqual('{:.3f}'.format(sum(seqrow)), '1.000')
        exp = [
            "TA\t0.125 0.176 0.143",
            "GG\t0.042 0.235 0.143",
            "AC\t0.167 0.059 0.143",
            "CT\t0.000 0.176 0.071",
            "AG\t0.042 0.118 0.071",
            "CA\t0.042 0.000 0.143",
            "AT\t0.042 0.059 0.071",
            "GA\t0.042 0.059 0.071",
            "AA\t0.083 0.000 0.071",
            "CC\t0.083 0.000 0.071",
            "CG\t0.125 0.000 0.000",
            "GT\t0.125 0.000 0.000",
            "TT\t0.042 0.059 0.000",
            "TC\t0.000 0.059 0.000",
            "TG\t0.042 0.000 0.000"
        ]
        self.assertEqual(freqs.format(), "\n".join(exp))

    def test_weightmodel_invalid_wtype(self):
        weights = {'A': 2, 'C': 2, 'G': 2, 'T': 2}
        with self.assertRaises(Exception) as context:
            weightmodel = word_vector.WeightModel(weights, 'nonexistent')        
        self.assertIn('unknown weight model', str(context.exception))

    def test_weighted_counts_pattern1(self):
        weights = {'A': 2, 'C': 2, 'G': 2, 'T': 2}
        weightmodel = word_vector.WeightModel(weights)
        cw = word_vector.CountsWeight(self.dna_records.length_list,
                                      self.pattern1, weightmodel)
        exp = ["A\t16 8 12",
               "G\t12 12 6",
               "C\t12 6 8",
               "T\t10 10 4"]
        self.assertEqual(cw.format(0), "\n".join(exp))
        for i in range(len(cw.data)):
            self.assertEqual(sum(cw.data[i]), cw.seq_lengths[i] * 2)

    def test_weighted_counts_pattern2(self):
        weights = {'A': 2, 'C': 2, 'G': 2, 'T': 2}
        weightmodel = word_vector.WeightModel(weights)
        cw = word_vector.CountsWeight(self.dna_records.length_list,
                                      self.pattern2, weightmodel)
        exp = [
            "TA\t12 12 8",
            "AC\t16 4 8",
            "GG\t4 16 8",
            "AG\t4 8 4",
            "CT\t0 12 4",
            "AA\t8 0 4",
            "AT\t4 4 4",
            "CA\t4 0 8",
            "CC\t8 0 4",
            "CG\t12 0 0",
            "GA\t4 4 4",
            "GT\t12 0 0",
            "TT\t4 4 0",
            "TC\t0 4 0",
            "TG\t4 0 0"
        ]
        self.assertEqual(cw.format(0), "\n".join(exp))

    def test_weighted_freqs_pattern1(self):
        weights = {'A': 2, 'C': 2, 'G': 2, 'T': 2}
        weightmodel = word_vector.WeightModel(weights)
        fw = word_vector.FreqsWeight(self.dna_records.length_list,
                                     self.pattern1, weightmodel)
        exp = [
            "A\t0.640 0.444 0.800",
            "G\t0.480 0.667 0.400",
            "C\t0.480 0.333 0.533",
            "T\t0.400 0.556 0.267"
        ]
        self.assertEqual(fw.format(), "\n".join(exp))

    def test_weighted_freqs_pattern1(self):
        weights = {'A': 2, 'C': 2, 'G': 2, 'T': 2}
        weightmodel = word_vector.WeightModel(weights)
        fw = word_vector.FreqsWeight(self.dna_records.length_list,
                                     self.pattern2, weightmodel)
        exp = [
            "TA\t0.500 0.706 0.571",
            "GG\t0.167 0.941 0.571",
            "AC\t0.667 0.235 0.571",
            "CT\t0.000 0.706 0.286",
            "AG\t0.167 0.471 0.286",
            "CA\t0.167 0.000 0.571",
            "AT\t0.167 0.235 0.286",
            "GA\t0.167 0.235 0.286",
            "AA\t0.333 0.000 0.286",
            "CC\t0.333 0.000 0.286",
            "CG\t0.500 0.000 0.000",
            "GT\t0.500 0.000 0.000",
            "TT\t0.167 0.235 0.000",
            "TC\t0.000 0.235 0.000",
            "TG\t0.167 0.000 0.000"
        ]
        self.assertEqual(fw.format(), "\n".join(exp))

    def test_equal_freqs_pattern2(self):
        # The result of this method is identical to that from decaf+py.
        p = word_pattern.create(self.dna_records.seq_list, 2, True)
        freq = word_vector.Freqs(self.dna_records.length_list, p)
        freqmodel = word_vector.EqualFreqs(alphabet_size=4)
        freqs_std = word_vector.FreqsStd(self.dna_records.length_list,
                                         p, freqmodel)
        exp = [
            "TA\t0.113 0.189 0.169",
            "AC\t0.150 0.063 0.169",
            "GG\t0.030 0.201 0.135",
            "CT\t0.000 0.189 0.084",
            "AG\t0.038 0.126 0.084",
            "CA\t0.038 0.000 0.169",
            "AT\t0.038 0.063 0.084",
            "GA\t0.038 0.063 0.084",
            "AA\t0.060 0.000 0.067",
            "CC\t0.060 0.000 0.067",
            "CG\t0.113 0.000 0.000",
            "GT\t0.113 0.000 0.000",
            "TT\t0.030 0.050 0.000",
            "TC\t0.000 0.063 0.000",
            "TG\t0.038 0.000 0.000"
        ]
        self.assertEqual(freqs_std.format(), "\n".join(exp))

    def test_equilibrium_freqs_pattern2(self):
        p = word_pattern.create(self.dna_records.seq_list, 2, True)
        dna_freqs = {'A': 0.24, 'C': 0.26, 'G': 0.23, 'T': 0.27}
        freqmodel = word_vector.EquilibriumFreqs(dna_freqs)
        freqs_std = word_vector.FreqsStd(self.dna_records.length_list,
                                         p, freqmodel)
        exp = [
            "TA\t0.111 0.186 0.166",
            "GG\t0.033 0.219 0.147",
            "AC\t0.151 0.063 0.169",
            "CT\t0.000 0.181 0.081",
            "AG\t0.040 0.132 0.089",
            "CA\t0.038 0.000 0.169",
            "GA\t0.040 0.066 0.089",
            "AT\t0.037 0.062 0.083",
            "AA\t0.062 0.000 0.070",
            "CC\t0.057 0.000 0.065",
            "CG\t0.115 0.000 0.000",
            "GT\t0.113 0.000 0.000",
            "TT\t0.028 0.046 0.000",
            "TC\t0.000 0.060 0.000",
            "TG\t0.038 0.000 0.000"
        ]
        self.assertEqual(freqs_std.format(), "\n".join(exp))

    def test_bools_pattern2(self):
        bools = word_vector.Bools(self.dna_records.length_list,
                                  self.pattern2)
        exp = [
            "AC\t1 1 1",
            "AG\t1 1 1",
            "AT\t1 1 1",
            "GA\t1 1 1",
            "GG\t1 1 1",
            "TA\t1 1 1",
            "AA\t1 0 1",
            "CA\t1 0 1",
            "CC\t1 0 1",
            "CT\t0 1 1",
            "TT\t1 1 0",
            "CG\t1 0 0",
            "GT\t1 0 0",
            "TC\t0 1 0",
            "TG\t1 0 0"
        ]
        self.assertEqual(bools.format(decimal_places=0), "\n".join(exp))

    def test_bools_pattern1(self):
        bools = word_vector.Bools(self.dna_records.length_list,
                                  self.pattern1)
        exp = [
            "A\t1 1 1",
            "C\t1 1 1",
            "G\t1 1 1",
            "T\t1 1 1"
        ]
        self.assertEqual(bools.format(decimal_places=0), "\n".join(exp))


if __name__ == '__main__':
    unittest.main()
