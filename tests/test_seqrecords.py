import os
import unittest
import uuid

from alfpy.utils import seqrecords


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


class TestSeqRecords(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def _validate_seqrecords(self, rec):
        self.assertEqual(rec.id_list, ID_LIST)
        self.assertEqual(rec.seq_list, SEQ_LIST)
        self.assertEqual(rec.length_list, [len(s) for s in SEQ_LIST])
        self.assertEqual(rec.count, len(SEQ_LIST))

    def test_SeqRecords_init(self):
        rec = seqrecords.SeqRecords(id_list=ID_LIST, seq_list=SEQ_LIST)
        self._validate_seqrecords(rec)

    def test_SeqRecords_add(self):
        rec = seqrecords.SeqRecords()
        for i in range(len(ID_LIST)):
            rec.add(ID_LIST[i], SEQ_LIST[i])
        self._validate_seqrecords(rec)

    def test_read_fasta(self):
        rec = seqrecords.read_fasta(FASTA_LIST)
        self._validate_seqrecords(rec)


if __name__ == '__main__':
    unittest.main()
