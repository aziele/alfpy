import os
import unittest

from alfpy import word_pattern

from . import utils


class Test(unittest.TestCase, utils.ModulesCommonTest):

    def __init__(self, *args, **kwargs):
        super(Test, self).__init__(*args, **kwargs)
        utils.ModulesCommonTest.set_test_data()

    def test_word_pattern_create_wordsize1_wordposFalse(self):
        p = word_pattern.create(self.dna_records.seq_list,
                                word_size=1,
                                wordpos=False)
        exp = [
            "18\t3\tA 0:8 1:4 2:6",
            "15\t3\tG 0:6 1:6 2:3",
            "13\t3\tC 0:6 1:3 2:4",
            "12\t3\tT 0:5 1:5 2:2"
        ]
        self.assertEqual(p.format(), "\n".join(exp))

    def test_word_pattern_create_wordsize1_wordposTrue(self):
        p = word_pattern.create(self.dna_records.seq_list,
                                word_size=1,
                                wordpos=True)
        exp = [
            '18\t3\tA 0 0 0 1 0 5 0 8 0 12 0 13 0 17 0 22 1 2 1 7 1 11 1 ' +
            '15 2 2 2 6 2 7 2 9 2 11 2 14',
            '15\t3\tG 0 3 0 11 0 15 0 20 0 23 0 24 1 3 1 4 1 5 1 6 1 16 1 ' +
            '17 2 3 2 4 2 5',
            '13\t3\tC 0 2 0 6 0 7 0 14 0 18 0 19 1 0 1 8 1 13 2 0 2 8 2 12 ' +
            '2 13',
            '12\t3\tT 0 4 0 9 0 10 0 16 0 21 1 1 1 9 1 10 1 12 1 14 2 ' +
            '1 2 10'
        ]
        self.assertEqual(p.format(), "\n".join(exp))

    def test_word_pattern_create_wordsize2_wordposFalse(self):
        p = word_pattern.create(self.dna_records.seq_list,
                                word_size=2,
                                wordpos=False)
        exp = ["8\t3\tTA 0:3 1:3 2:2",
               "7\t3\tAC 0:4 1:1 2:2",
               "7\t3\tGG 0:1 1:4 2:2",
               "4\t3\tAG 0:1 1:2 2:1",
               "4\t2\tCT 1:3 2:1",
               "3\t3\tAT 0:1 1:1 2:1",
               "3\t3\tGA 0:1 1:1 2:1",
               "3\t2\tAA 0:2 2:1",
               "3\t2\tCA 0:1 2:2",
               "3\t2\tCC 0:2 2:1",
               "3\t1\tCG 0:3",
               "3\t1\tGT 0:3",
               "2\t2\tTT 0:1 1:1",
               "1\t1\tTC 1:1",
               "1\t1\tTG 0:1",
               ]
        self.assertEqual(p.format(), "\n".join(exp))

    def test_input_output_file_pattern(self):

        for wordpos in [True, False]:
            p1 = word_pattern.create(self.dna_records.seq_list,
                                     word_size=1,
                                     wordpos=wordpos)
            oh = open(utils.get_test_data('pattern.txt'), 'w')
            oh.write(p1.format())
            oh.close()
            fh = open(utils.get_test_data('pattern.txt'))
            p2 = word_pattern.read(fh)
            fh.close()
            self.assertEqual(p1.format(), p2.format())
        os.remove(utils.get_test_data('pattern.txt'))

    def test_reduce_alphabet_wordsize1(self):
        alphabet_dict = {'A': 'R', 'C': 'Y', 'T': 'Y', 'G': 'R'}
        p = word_pattern.create(self.dna_records.seq_list,
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
        p = word_pattern.create(self.dna_records.seq_list,
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
        p = word_pattern.create(self.dna_records.seq_list,
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

    def test_create_from_fasta_dna_word_size2(self):
        fh = open(self.dna_filename)
        p = word_pattern.create_from_fasta(fh, word_size=2)
        fh.close()
        md5 = utils.calc_md5(p.format())
        self.assertEqual(md5, '5be104951119e2df4528100b7fc672f4')

    def test_create_from_bigfasta_dna_word_size2(self):
        p = word_pattern.create_from_bigfasta(self.dna_filename, 2)
        md5 = utils.calc_md5(p.format())
        self.assertEqual(md5, '5be104951119e2df4528100b7fc672f4')

    def test_create_from_fasta_pep_word_size1(self):
        fh = open(self.pep_filename)
        p = word_pattern.create_from_fasta(fh, word_size=1)
        fh.close()
        md5 = utils.calc_md5(p.format())
        self.assertEqual(md5, '2d4dd98798cb6320975f6919fe43b777')

    def test_create_from_bigfasta_pep_word_size1(self):
        p = word_pattern.create_from_bigfasta(self.pep_filename, 1)
        md5 = utils.calc_md5(p.format())
        self.assertEqual(md5, '2d4dd98798cb6320975f6919fe43b777')


if __name__ == '__main__':
    unittest.main()
