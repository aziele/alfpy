import unittest

from alfpy.utils import seqrecords

from . import utils


class SeqRecordsTest(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(SeqRecordsTest, self).__init__(*args, **kwargs)
        self.ID_LIST = ['seq1', 'seq2', 'seq3', 'seq4']
        self.DESC_LIST = ['seq1 desc', 'seq2 desc', 'seq3 desc', '']
        self.SEQ_LIST = [
            'MEVVIRSANFTDNAKIIIVQLNASVEINCTRPNNYTRKGIRIGPGRAVYAAEEIIGDNTLKQVVTKLRE',
            'MVIRSANFTDNAKIIIVQLNASVEINCTRPNNNTRKGIRIGPGRAVYAAEEIIGDIRRAHCNIS',
            'MFTDNAKIIIVQLNASVEINCTRPNNNTRKGIHIGPGRAFYATGEIIGDIRQAHCNISGAKW',
            'MFTDNAKIIIVQLNASVEINCTRPNNNTR'
        ]

    def _validate_seqrecords(self, rec):
        self.assertEqual(rec.id_list, self.ID_LIST)
        self.assertEqual(rec.seq_list, self.SEQ_LIST)
        self.assertEqual(rec.length_list, [len(s) for s in self.SEQ_LIST])
        self.assertEqual(rec.count, len(self.SEQ_LIST))

    def test_SeqRecords_init(self):
        rec = seqrecords.SeqRecords(
            id_list=self.ID_LIST, seq_list=self.SEQ_LIST)
        self._validate_seqrecords(rec)

    def test_SeqRecords_add(self):
        rec = seqrecords.SeqRecords()
        for i in range(len(self.ID_LIST)):
            rec.add(self.ID_LIST[i], self.SEQ_LIST[i])
        self._validate_seqrecords(rec)

    def test_read_fasta(self):
        fh = open(utils.get_test_data('pep.fa'))
        rec = seqrecords.read_fasta(fh)
        fh.close()
        self._validate_seqrecords(rec)


if __name__ == '__main__':
    unittest.main()
