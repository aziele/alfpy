
import os
import unittest
import uuid

from alfpy.utils import seqrecords
from alfpy import word_pattern
from alfpy import word_vector


ID_LIST = ['seq1', 'seq2', 'seq3', 'seq4']
SEQ_LIST = [
    'AACGTACCATTGAACGTACCGTAGG',
    'ctaggggacttatctagg',
    'CTAGGGAACATACCA'
]


class TestWordVector(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestWordVector, self).__init__(*args, **kwargs)
        self.seq_records = seqrecords.SeqRecords(ID_LIST, SEQ_LIST)
        self.pattern1 = word_pattern.create(self.seq_records.seq_list, 1)
        self.pattern2 = word_pattern.create(self.seq_records.seq_list, 2)

    def test_counts_pattern1(self):
        counts = word_vector.Counts(self.seq_records.length_list, 
                                    self.pattern1)
        pat_list =  ['A', 'C', 'T', 'G']
        data = [[8, 6, 5, 6],
                [4, 3, 5, 6],
                [6, 4, 2, 3]]
        self.assertEqual(counts.pat_list, pat_list)  
        self.assertEqual(counts.data.tolist(), list(data))
        # Check if counts sum to sequences length
        for i in range(len(data)):
            self.assertEqual(sum(data[i]), counts.seq_lengths[i])

    def test_counts_pattern2(self):
        counts = word_vector.Counts(self.seq_records.length_list, 
                                    self.pattern2)
        pat_list = ['AA', 'AC', 'GT', 'AG', 'CC', 'CA', 'CG', 'TC', 'GG', 
                    'AT', 'GA', 'TG', 'CT', 'TT', 'TA']
        data = [[2, 4, 3, 1, 2, 1, 3, 0, 1, 1, 1, 1, 0, 1, 3], 
                [0, 1, 0, 2, 0, 0, 0, 1, 4, 1, 1, 0, 3, 1, 3], 
                [1, 2, 0, 1, 1, 2, 0, 0, 2, 1, 1, 0, 1, 0, 2]]
        self.assertEqual(counts.pat_list, pat_list)  
        self.assertEqual(counts.data.tolist(), list(data))
        # Check if counts sum to sequences length
        for i in range(len(data)):
            self.assertEqual(sum(data[i]), counts.seq_lengths[i]-1)


    def test_freqs_pattern1(self):
        freqs = word_vector.Freqs(self.seq_records.length_list, 
                                 self.pattern1)
        pat_list =  ['A', 'C', 'T', 'G']
        # Check whether freqs in a given sequence sum to 1.
        for seqrow in freqs.data:
            self.assertEqual('{:.3f}'.format(sum(seqrow)), '1.000')
        data = ['0.320\t0.240\t0.200\t0.240',
                '0.222\t0.167\t0.278\t0.333',      
                '0.400\t0.267\t0.133\t0.200']
        self.assertEqual(freqs.format(), "\n".join(data))

    def test_freqs_pattern2(self):
        freqs = word_vector.Freqs(self.seq_records.length_list, 
                                 self.pattern2)
        pat_list = ['AA', 'AC', 'GT', 'AG', 'CC', 'CA', 'CG', 'TC', 'GG', 
                    'AT', 'GA', 'TG', 'CT', 'TT', 'TA']
        # Check whether freqs in a given sequence sum to 1.
        for seqrow in freqs.data:
            self.assertEqual('{:.3f}'.format(sum(seqrow)), '1.000')
        data = ['0.083\t0.167\t0.125\t0.042\t0.083\t0.042\t0.125\t0.000' +
                '\t0.042\t0.042\t0.042\t0.042\t0.000\t0.042\t0.125',
                '0.000\t0.059\t0.000\t0.118\t0.000\t0.000\t0.000\t0.059' +
                '\t0.235\t0.059\t0.059\t0.000\t0.176\t0.059\t0.176',
                '0.071\t0.143\t0.000\t0.071\t0.071\t0.143\t0.000\t0.000' +
                '\t0.143\t0.071\t0.071\t0.000\t0.071\t0.000\t0.143']
        self.assertEqual(freqs.format(), "\n".join(data)) 


    def test_weighted_counts_pattern1(self):
        weights = {'A': 2, 'C': 2, 'G': 2, 'T': 2}
        weightmodel = word_vector.WeightModel(weights)
        cw = word_vector.CountsWeight(self.seq_records.length_list, 
                                          self.pattern1, weightmodel)
        pat_list = ['A', 'C', 'T', 'G']
        data = [[16, 12, 10, 12],
                [8, 6, 10, 12],
                [12, 8, 4, 6]]
        self.assertEqual(cw.pat_list, pat_list)  
        self.assertEqual(cw.data.tolist(), list(data))

    def test_weighted_counts_pattern2(self):
        weights = {'A': 2, 'C': 2, 'G': 2, 'T': 2}
        weightmodel = word_vector.WeightModel(weights)
        cw = word_vector.CountsWeight(self.seq_records.length_list, 
                                          self.pattern2, weightmodel)
        pat_list = ['AA', 'AC', 'GT', 'AG', 'CC', 'CA', 'CG', 'TC', 'GG', 
                    'AT', 'GA', 'TG', 'CT', 'TT', 'TA']
        data = [[8, 16, 12, 4, 8, 4, 12, 0, 4, 4, 4, 4, 0, 4, 12],
                [0, 4, 0, 8, 0, 0, 0, 4, 16, 4, 4, 0, 12, 4, 12],
                [4, 8, 0, 4, 4, 8, 0, 0, 8, 4, 4, 0, 4, 0, 8]]
        self.assertEqual(cw.pat_list, pat_list)  
        self.assertEqual(cw.data.tolist(), data)


    def test_weighted_freqs_pattern1(self):
        weights = {'A': 2, 'C': 2, 'G': 2, 'T': 2}
        weightmodel = word_vector.WeightModel(weights)
        fw = word_vector.FreqsWeight(self.seq_records.length_list, 
                                          self.pattern1, weightmodel)
        pat_list = ['A', 'C', 'T', 'G']
        data = ['0.640\t0.480\t0.400\t0.480',
                '0.444\t0.333\t0.556\t0.667',
                '0.800\t0.533\t0.267\t0.400']
        self.assertEqual(fw.format(), "\n".join(data))

    def test_weighted_freqs_pattern1(self):
        weights = {'A': 2, 'C': 2, 'G': 2, 'T': 2}
        weightmodel = word_vector.WeightModel(weights)
        fw = word_vector.FreqsWeight(self.seq_records.length_list, 
                                          self.pattern2, weightmodel)
        pat_list = ['AA', 'AC', 'GT', 'AG', 'CC', 'CA', 'CG', 'TC', 'GG', 
                    'AT', 'GA', 'TG', 'CT', 'TT', 'TA']
        data = ['0.333\t0.667\t0.500\t0.167\t0.333\t0.167\t0.500\t0.000' + 
                '\t0.167\t0.167\t0.167\t0.167\t0.000\t0.167\t0.500',
                '0.000\t0.235\t0.000\t0.471\t0.000\t0.000\t0.000\t0.235' + 
                '\t0.941\t0.235\t0.235\t0.000\t0.706\t0.235\t0.706',
                '0.286\t0.571\t0.000\t0.286\t0.286\t0.571\t0.000\t0.000' +
                '\t0.571\t0.286\t0.286\t0.000\t0.286\t0.000\t0.571']
        self.assertEqual(fw.format(), "\n".join(data))


    def test_weighted_freqs_pattern1(self):
        weights = {'A': 2, 'C': 2, 'G': 2, 'T': 2}
        weightmodel = word_vector.WeightModel(weights)
        fw = word_vector.FreqsWeight(self.seq_records.length_list, 
                                          self.pattern2, weightmodel)
        pat_list = ['AA', 'AC', 'GT', 'AG', 'CC', 'CA', 'CG', 'TC', 'GG', 
                    'AT', 'GA', 'TG', 'CT', 'TT', 'TA']
        data = ['0.333\t0.667\t0.500\t0.167\t0.333\t0.167\t0.500\t0.000' + 
                '\t0.167\t0.167\t0.167\t0.167\t0.000\t0.167\t0.500',
                '0.000\t0.235\t0.000\t0.471\t0.000\t0.000\t0.000\t0.235' + 
                '\t0.941\t0.235\t0.235\t0.000\t0.706\t0.235\t0.706',
                '0.286\t0.571\t0.000\t0.286\t0.286\t0.571\t0.000\t0.000' +
                '\t0.571\t0.286\t0.286\t0.000\t0.286\t0.000\t0.571']
        self.assertEqual(fw.format(), "\n".join(data))

    def test_equilibrium_freqs_pattern2(self):
        p = word_pattern.create(self.seq_records.seq_list, 2, True)
        freq = word_vector.Freqs(self.seq_records.length_list, p)
        freqmodel = word_vector.EqualFreqs(alphabet_size=4)
        freqs_std = word_vector.FreqsStd(self.seq_records.length_list, 
                                         p, freqmodel)
        # Freqs_std agrees with decaf+py
        pat_list = ['AA', 'AC', 'GT', 'AG', 'CC', 'CA', 'CG', 'TC', 'GG', 
                    'AT', 'GA', 'TG', 'CT', 'TT', 'TA']
        data = ['0.060\t0.150\t0.113\t0.038\t0.060\t0.038\t0.113\t0.000' + 
                '\t0.030\t0.038\t0.038\t0.038\t0.000\t0.030\t0.113',
                '0.000\t0.063\t0.000\t0.126\t0.000\t0.000\t0.000\t0.063' +
                '\t0.201\t0.063\t0.063\t0.000\t0.189\t0.050\t0.189',
                '0.067\t0.169\t0.000\t0.084\t0.067\t0.169\t0.000\t0.000' +
                '\t0.135\t0.084\t0.084\t0.000\t0.084\t0.000\t0.169']
        self.assertEqual(freqs_std.format(), "\n".join(data))



if __name__ == '__main__':
    unittest.main()
