import os
import unittest
import uuid

from alfpy.utils import fasta


# Input data for tests.
ID_LIST = ['seq1', 'seq2', 'seq3', 'seq4']
DESC_LIST = ['seq1 desc', 'seq2 desc', 'seq3 desc', '']
SEQ_LIST = [
    'MEVVIRSANFTDNAKIIIVQLNASVEINCTRPNNYTRKGIRIGPGRAVYAAEEIIGDNTLKQVVTKLRE',
    'MVIRSANFTDNAKIIIVQLNASVEINCTRPNNNTRKGIRIGPGRAVYAAEEIIGDIRRAHCNIS',
    'MFTDNAKIIIVQLNASVEINCTRPNNNTRKGIHIGPGRAFYATGEIIGDIRQAHCNISGAKW',
    'MFTDNAKIIIVQLNASVEINCTRPNNNTR'
]
FASTA_LIST = ['>seq1 seq1 desc\n',
              'MEVVIRSANFTDNAKIIIVQLNASVEINCTR\n',
              'PNNYTRKGIRIGPGRAVYAAEEIIGDNTLKQ\n',
              'VVTKLRE\n',
              '>seq2 seq2 desc\n',
              'MVIRSANFTDNAKIIIVQLNASVEINCTRPN\n',
              'NNTRKGIRIGPGRAVYAAEEIIGDIRRAHCN\n',
              'IS\n',
              '>seq3 seq3 desc\n',
              'MFTDNAKIIIVQLNASVEINCTRPNNNTRKG\n',
              'IHIGPGRAFYATGEIIGDIRQAHCNISGAKW\n',
              '>seq4\n',
              'MFTDNAKIIIVQLNASVEINCTRPNNNTR\n'
              ]
#


class TestFasta(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def _validate_FastaRecord_init(self, fasta_record, seqidx):
        self.assertIsInstance(fasta_record, fasta.FastaRecord)
        self.assertEqual(fasta_record.seq, SEQ_LIST[seqidx])
        self.assertEqual(fasta_record.id, ID_LIST[seqidx])
        self.assertEqual(fasta_record.description, DESC_LIST[seqidx])
        self.assertEqual(len(fasta_record), len(SEQ_LIST[seqidx]))

    def test_single_FastaRecord_init(self):
        r = fasta.FastaRecord(SEQ_LIST[0], ID_LIST[0], DESC_LIST[0])
        self._validate_FastaRecord_init(r, seqidx=0)

    def test_multiple_FastaRecord_init(self):
        for i in range(len(ID_LIST)):
            r = fasta.FastaRecord(SEQ_LIST[i], ID_LIST[i], DESC_LIST[i])
            self._validate_FastaRecord_init(r, seqidx=i)

    def test_read_fasta(self):
        r = fasta.read(FASTA_LIST)
        self._validate_FastaRecord_init(r, seqidx=0)

    def test_parse_fasta(self):
        for i, rec in enumerate(fasta.parse(FASTA_LIST)):
            self._validate_FastaRecord_init(rec, seqidx=i)

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
        filename = '{}'.format(uuid.uuid4().hex)
        oh = open(filename, 'w')
        l1 = []
        for seq_record in fasta.parse(FASTA_LIST):
            l1.append(seq_record.format())
            oh.write(seq_record.format())
            oh.write('\n')
        oh.close()
        fh = open(filename)
        l2 = [seq_record.format() for seq_record in fasta.parse(fh)]
        fh.close()
        os.remove(filename)
        self.assertEqual(l1, l2)


if __name__ == '__main__':
    unittest.main()
