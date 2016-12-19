import unittest

from . import utils


class ScriptTest(unittest.TestCase, utils.ScriptsCommonTest):

    def __init__(self, *args, **kwargs):
        super(ScriptTest, self).__init__(*args, **kwargs)
        utils.ScriptsCommonTest.set_test_data()
        self.script_name = 'calc_fcgr.py'

    def test_arg_word_size_2_when_no_fasta(self):
        args = ['--word_size', '2']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--fasta/-f', out)

    def test_arg_fasta_when_no_word_size(self):
        args = ['--fasta', self.filename_dna]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--word_size/-w', out)

    def test_arg_word_size_too_small(self):
        args = ['--fasta', self.filename_dna, '--word_size', '0']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--word_size must be >= 1', out)

    def test_output_word_size_1(self):
        args = ['--fasta', self.filename_dna, '--word_size', '1']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, 'bee51f3214f06f4e4265aa05bf9d6a7e')

    def test_output_word_size_2(self):
        args = ['--fasta', self.filename_dna, '--word_size', '2']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '7175a91fb9fc31661ce07aea28743605')


if __name__ == '__main__':
    unittest.main()
