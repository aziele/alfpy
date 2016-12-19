import unittest

from . import utils


class ScriptTest(unittest.TestCase, utils.ScriptsWordCommonTest):

    def __init__(self, *args, **kwargs):
        super(ScriptTest, self).__init__(*args, **kwargs)
        utils.ScriptsWordCommonTest.set_test_data()
        self.script_name = 'calc_word_bool.py'

    def test_arg_word_size_when_no_fasta(self):
        args = ['--word_size', '2']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--fasta/-f', out)

    def test_arg_word_pattern_when_no_fasta(self):
        args = ['--word_pattern', self.filename_pep_2mer]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--fasta/-f', out)

    def test_arg_fasta_when_no_wordsize_or_wordpattern(self):
        args = ['--fasta', self.filename_pep]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Specify either: --word_size or --word', out)

    def test_arg_fasta_when_no_wordsize_or_wordpattern(self):
        args = ['--fasta', self.filename_pep]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Specify either: --word_size or --word', out)

    def test_arg_word_size_too_small(self):
        args = ['--fasta', self.filename_pep, '--word_size', '-1']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Word size must be >= 1.', out)

    def test_output_word_size1(self):
        args = ['--fasta', self.filename_pep, '--word_size', '1']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '4caed60c7590f45e9a6de19482839e9c')


if __name__ == '__main__':
    unittest.main()
