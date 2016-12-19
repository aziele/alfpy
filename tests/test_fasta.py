import os
import unittest

from alfpy.utils import fasta

from . import utils


class FastaTest(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(FastaTest, self).__init__(*args, **kwargs)
        self.ID_LIST = ['seq1', 'seq2', 'seq3', 'seq4']
        self.DESC_LIST = ['seq1 desc', 'seq2 desc', 'seq3 desc', '']
        self.SEQ_LIST = [
         'MEVVIRSANFTDNAKIIIVQLNASVEINCTRPNNYTRKGIRIGPGRAVYAAEEIIGDNTLKQVVTKLRE',
         'MVIRSANFTDNAKIIIVQLNASVEINCTRPNNNTRKGIRIGPGRAVYAAEEIIGDIRRAHCNIS',
         'MFTDNAKIIIVQLNASVEINCTRPNNNTRKGIHIGPGRAFYATGEIIGDIRQAHCNISGAKW',
         'MFTDNAKIIIVQLNASVEINCTRPNNNTR'
        ]

    def _validate_FastaRecord_init(self, fasta_record, seqidx):
        self.assertIsInstance(fasta_record, fasta.FastaRecord)
        self.assertEqual(fasta_record.seq, self.SEQ_LIST[seqidx])
        self.assertEqual(fasta_record.id, self.ID_LIST[seqidx])
        self.assertEqual(fasta_record.description, self.DESC_LIST[seqidx])
        self.assertEqual(len(fasta_record), len(self.SEQ_LIST[seqidx]))

    def test_single_FastaRecord_init(self):
        r = fasta.FastaRecord(self.SEQ_LIST[0],
                              self.ID_LIST[0],
                              self.DESC_LIST[0])
        self._validate_FastaRecord_init(r, seqidx=0)

    def test_multiple_FastaRecord_init(self):
        for i in range(len(self.ID_LIST)):
            r = fasta.FastaRecord(self.SEQ_LIST[i],
                                  self.ID_LIST[i],
                                  self.DESC_LIST[i])
            self._validate_FastaRecord_init(r, seqidx=i)

    def test_read_fasta(self):
        fh = open(utils.get_test_data('pep.fa'))
        r = fasta.read(fh)
        fh.close()
        self._validate_FastaRecord_init(r, seqidx=0)

    def test_parse_fasta(self):
        fh = open(utils.get_test_data('pep.fa'))
        for i, rec in enumerate(fasta.parse(fh)):
            self._validate_FastaRecord_init(rec, seqidx=i)
        fh.close()

    def test_parse_fasta_missing_sequences(self):
        ids = ['seq1', 'seq2']
        seqs = ['ATGC', '']
        l = ['>{}\n'.format(ids[0]),
             '{}\n\n\n'.format(seqs[0]),
             '>{}\n'.format(ids[1]),
             '{}\n'.format(seqs[1])
             ]
        for i, fasta_record in enumerate(fasta.parse(l)):
            self.assertIsInstance(fasta_record, fasta.FastaRecord)
            self.assertEqual(fasta_record.seq, seqs[i])

    def test_fasta_format(self, wrap=70):
        l = ['>seq1 seq1 desc\n',
             'A' * wrap + '\n',
             'B' * wrap]
        r = fasta.read(l)
        self.assertEqual(''.join(l), r.format(wrap=wrap))

    def test_input_output_file_fasta(self):
        filename = 'temp.fa'
        oh = open(utils.get_test_data(filename), 'w')
        l1 = []
        fh = open(utils.get_test_data('pep.fa'))
        for seq_record in fasta.parse(fh):
            l1.append(seq_record.format())
            oh.write(seq_record.format())
            oh.write('\n')
        fh.close()
        oh.close()
        fh = open(utils.get_test_data(filename))
        l2 = [seq_record.format() for seq_record in fasta.parse(fh)]
        fh.close()
        os.remove(utils.get_test_data(filename))
        self.assertEqual(l1, l2)


if __name__ == '__main__':
    unittest.main()
