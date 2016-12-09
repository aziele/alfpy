import os
import unittest
import uuid

from alfpy.utils import seqrecords
from alfpy import word_pattern


FASTA = ['>seq1 seq1 desc\n',
         'AACGTACCATTGAACGTACCGTAGG\n',
         '>seq2 seq2 desc\n',
         'CTAGGGGACTTATCTAGG\n',
         '>seq3 seq3 desc\n',
         'CTAGGGAACATACCA\n'
         ]
PATTERN1 = {
    'pat_list': ['A', 'C', 'T', 'G'],
    'occr_list': [
        {0: 8, 1: 4, 2: 6},
        {0: 6, 1: 3, 2: 4},
        {0: 5, 1: 5, 2: 2},
        {0: 6, 1: 6, 2: 3}
    ],
    'pos_list': [
        {0: [0, 1, 5, 8, 12, 13, 17, 22],
         1: [2, 7, 11, 15],
         2: [2, 6, 7, 9, 11, 14]},
        {0: [2, 6, 7, 14, 18, 19],
            1: [0, 8, 13],
            2: [0, 8, 12, 13]},
        {0: [4, 9, 10, 16, 21],
            1: [1, 9, 10, 12, 14],
            2: [1, 10]},
        {0: [3, 11, 15, 20, 23, 24],
            1: [3, 4, 5, 6, 16, 17],
            2: [3, 4, 5]}],
}
PATTERN2 = {
    'pat_list': [
        'AA', 'AC', 'GT', 'AG', 'CC', 'CA', 'CG', 'TC',
        'GG', 'AT', 'GA', 'TG', 'CT', 'TT', 'TA'
    ],
    'occr_list': [
        {0: 2, 2: 1},
        {0: 4, 1: 1, 2: 2},
        {0: 3},
        {0: 1, 1: 2, 2: 1},
        {0: 2, 2: 1},
        {0: 1, 2: 2},
        {0: 3},
        {1: 1},
        {0: 1, 1: 4, 2: 2},
        {0: 1, 1: 1, 2: 1},
        {0: 1, 1: 1, 2: 1},
        {0: 1},
        {1: 3, 2: 1},
        {0: 1, 1: 1},
        {0: 3, 1: 3, 2: 2}
    ]
}


class TestPattern(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestPattern, self).__init__(*args, **kwargs)
        self.seq_records = seqrecords.read_fasta(FASTA)

    def test_word_pattern_create_wordsize1_wordposFalse(self):
        p = word_pattern.create(self.seq_records.seq_list,
                                word_size=1,
                                wordpos=False)
        self.assertEqual(p.pat_list, PATTERN1['pat_list'])
        self.assertEqual(p.occr_list, PATTERN1['occr_list'])
        self.assertEqual(p.pos_list, [])

    def test_word_pattern_create_wordsize1_wordposTrue(self):
        p = word_pattern.create(self.seq_records.seq_list,
                                word_size=1,
                                wordpos=True)
        self.assertEqual(p.pat_list, PATTERN1['pat_list'])
        self.assertEqual(p.occr_list, PATTERN1['occr_list'])
        self.assertEqual(p.pos_list, PATTERN1['pos_list'])

    def test_word_pattern_create_wordsize2_wordposFalse(self):
        p = word_pattern.create(self.seq_records.seq_list,
                                word_size=2,
                                wordpos=False)
        self.assertEqual(p.pat_list, PATTERN2['pat_list'])
        self.assertEqual(p.occr_list, PATTERN2['occr_list'])

    def test_pattern_format_wordsize1_wordposFalse(self):
        p = word_pattern.create(self.seq_records.seq_list,
                                word_size=1,
                                wordpos=False)
        expected_format = [
            '18\t3\tA 0:8 1:4 2:6',
            '15\t3\tG 0:6 1:6 2:3',
            '13\t3\tC 0:6 1:3 2:4',
            '12\t3\tT 0:5 1:5 2:2'
        ]
        self.assertEqual(p.format(), "\n".join(expected_format))

    def test_pattern_format_wordsize1_wordposTrue(self):
        p = word_pattern.create(self.seq_records.seq_list,
                                word_size=1,
                                wordpos=True)
        expected_format = [
            '18\t3\tA 0 0 0 1 0 5 0 8 0 12 0 13 0 17 0 22 1 2 1 7 1 11 1 ' +
            '15 2 2 2 6 2 7 2 9 2 11 2 14',
            '15\t3\tG 0 3 0 11 0 15 0 20 0 23 0 24 1 3 1 4 1 5 1 6 1 16 1 ' +
            '17 2 3 2 4 2 5',
            '13\t3\tC 0 2 0 6 0 7 0 14 0 18 0 19 1 0 1 8 1 13 2 0 2 8 2 12 ' +
            '2 13',
            '12\t3\tT 0 4 0 9 0 10 0 16 0 21 1 1 1 9 1 10 1 12 1 14 2 ' +
            '1 2 10'
        ]
        self.assertEqual(p.format(), "\n".join(expected_format))

    def test_input_output_file_pattern(self):
        filename = '{}.pattern'.format(uuid.uuid4().hex)
        for wordpos in [True, False]:
            p1 = word_pattern.create(self.seq_records.seq_list,
                                     word_size=1,
                                     wordpos=wordpos)
            oh = open(filename, 'w')
            oh.write(p1.format())
            oh.close()
            fh = open(filename)
            p2 = word_pattern.read(fh)
            fh.close()
            self.assertEqual(p1.format(), p2.format())
        os.remove(filename)

    def test_reduce_alphabet_wordsize1(self):
        alphabet_dict = {'A': 'R', 'C': 'Y', 'T': 'Y', 'G': 'R'}
        p = word_pattern.create(self.seq_records.seq_list,
                                word_size=1,
                                wordpos=False)

        expected_format = [
            '33\t3\tR 0:14 1:10 2:9',
            '25\t3\tY 0:11 1:8 2:6'
        ]
        p = p.reduce_alphabet(alphabet_dict)
        self.assertEqual(p.format(), "\n".join(expected_format))

    def test_reduce_alphabet_wordsize2(self):
        alphabet_dict = {'A': 'R', 'C': 'Y', 'T': 'Y', 'G': 'R'}
        p = word_pattern.create(self.seq_records.seq_list,
                                word_size=2,
                                wordpos=False)
        expected_format = [
            '17\t3\tRR 0:5 1:7 2:5',
            '15\t3\tYR 0:8 1:3 2:4',
            '13\t3\tRY 0:8 1:2 2:3',
            '10\t3\tYY 0:3 1:5 2:2'
        ]
        p = p.reduce_alphabet(alphabet_dict)
        self.assertEqual(p.format(), "\n".join(expected_format))

    def test_reduce_alphabet_wordsize1(self):
        p = word_pattern.create(self.seq_records.seq_list,
                                word_size=2,
                                wordpos=False)
        p1 = p.merge_revcomp()
        pat_list = ['AA', 'AC', 'AG', 'CC', 'CA', 'CG', 'AT', 'GA', 'TA']
        occr_list = [
            {0: 3, 1: 1, 2: 1},
            {0: 7, 1: 1, 2: 2},
            {0: 1, 1: 5, 2: 2},
            {0: 3, 1: 4, 2: 3},
            {0: 2, 2: 2},
            {0: 3},
            {0: 1, 1: 1, 2: 1},
            {0: 1, 1: 2, 2: 1},
            {0: 3, 1: 3, 2: 2}
        ]
        p2 = word_pattern.Pattern(pat_list, occr_list, [])
        self.assertEqual(p1.format(), p2.format())


if __name__ == '__main__':
    unittest.main()
